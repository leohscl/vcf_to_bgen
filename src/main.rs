use bgen_reader::bgen::bgen_stream::write_samples;
use bgen_reader::bgen::header::{Header, HeaderFlags};
use bgen_reader::bgen::variant_data::{DataBlock, VariantData};
use color_eyre::eyre::Report;
use color_eyre::Result;
use flate2::read::MultiGzDecoder;
use indicatif::ProgressBar;
use nom::bytes::complete::{is_not, tag, take, take_until, take_while1};
use nom::character::complete::{alpha0, alphanumeric0, char, tab};
use nom::multi::{many0, separated_list0};
use nom::sequence::{preceded, terminated};
use nom::{IResult, InputIter};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::time::Duration;

fn main() -> Result<()> {
    // let mut reader = BufReader::new(MultiGzDecoder::new(File::open(
    //     "/pasteur/zeus/projets/p02/GGS_WKD/HOME_LEO/vcf_to_bgen/data/small.vcf.gz"
    // )?));
    let mut reader = BufReader::new(MultiGzDecoder::new(File::open(
        "data/ukb24304_c1_b62_v1.vcf.gz",
    )?));
    // First pass to get the number of variants
    let mut variant_num = 0;
    let mut line = String::new();
    let bar = ProgressBar::new_spinner();
    bar.enable_steady_tick(Duration::from_millis(100));
    println!("Counting variants...  ");
    loop {
        let num_bytes = reader.read_line(&mut line)?;
        if num_bytes == 0 {
            break;
        }
        if !line.starts_with('#') {
            variant_num += 1;
        }
        line.clear();
    }
    bar.finish();
    println!("Done");
    drop(reader);

    // read file again
    let mut reader = BufReader::new(MultiGzDecoder::new(File::open(
        "data/ukb24304_c1_b62_v1.vcf.gz",
    )?));

    // Skip header, parse column/sample line
    let samples_result = loop {
        reader.read_line(&mut line)?;
        if line.starts_with("##") {
            line.clear();
            continue;
        } else if line.starts_with('#') {
            break parse_samples(&line);
        }
    };
    let samples_str = samples_result
        .map_err(|_err| Report::msg("Unable to parse samples line in header"))?
        .1;
    let samples: Vec<String> = samples_str.into_iter().map(|s| s.to_string()).collect();
    let number_individuals = samples.len() as u32;
    line.clear();

    // write bgen file
    let mut bgen_writer = BufWriter::new(File::create("data/ukb24304_c1_b62_v1.bgen")?);
    // compute length of sample block
    let len_sample_block =
        8u32 + number_individuals * 2 + samples.iter().map(|s| s.len() as u32).sum::<u32>();

    // compute length of header
    let header_size = 20u32;

    // compute offset to start of data
    let start_data_offset = header_size + len_sample_block;

    // create bgen header
    let header_flags = HeaderFlags {
        compressed_snp_blocks: true,
        layout_id: 2,
        sample_id_present: true,
    };
    let header = Header {
        start_data_offset,
        header_size,
        variant_num,
        variant_count: 0,
        sample_num: number_individuals,
        header_flags,
    };
    // write header
    header.write_header(&mut bgen_writer)?;

    // write samples
    write_samples(&samples, &mut bgen_writer, len_sample_block)?;

    // write variant blocks
    println!("Converting variants to bgen format");
    let bar = ProgressBar::new(variant_num as u64);
    for _variant_i in 0..variant_num {
        reader.read_line(&mut line)?;
        let variant_data = parse_genotype_line(&line, number_individuals)?;
        variant_data.write_self(&mut bgen_writer, 2)?;
        bar.inc(1);
        line.clear();
    }
    bar.finish();

    Ok(())
}

fn parse_samples(input: &str) -> IResult<&str, Vec<&str>> {
    preceded(
        preceded(tag("#"), many0(preceded(alpha0, tab))),
        separated_list0(tab, alphanumeric0),
    )(input)
}

fn parse_one_field(input: &str) -> Result<(&str, &str)> {
    terminated(is_not("\t"), char('\t'))(input).map_err(
        |_err: nom::Err<nom::error::Error<&str>>| {
            Report::msg("Unable to parse info fields of genotype")
        },
    )
}

fn parse_genotype_line(input: &str, number_individuals: u32) -> Result<VariantData> {
    let (remaining_input, chr) = parse_one_field(input)?;
    let (remaining_input, pos) = parse_one_field(remaining_input)?;
    let (remaining_input, variant_id) = parse_one_field(remaining_input)?;
    let (remaining_input, a1) = parse_one_field(remaining_input)?;
    let (remaining_input, a2) = parse_one_field(remaining_input)?;
    let genos_string = parse_genotype_field(remaining_input)
        .map_err(|_err| Report::msg("Unable to parse genotype information"))?
        .1;
    let variant_id_fmt = format_id_with_alleles(variant_id, a1, a2)?;
    // dbg!(&genos_string);
    let ploidy_missingness = genos_string
        .iter()
        .map(|geno_s| {
            // 2 if for ploidy
            if geno_s.contains('.') {
                (1u8 << 7) + 2
            } else {
                2u8
            }
        })
        .collect();

    let probabilities = genos_string
        .into_iter()
        .flat_map(|geno_s| {
            let mut geno_iter = geno_s.iter_elements().filter_map(|c| c.to_digit(10));
            let left_strand = geno_iter.next().unwrap_or(0);
            let right_strand = geno_iter.nth(1).unwrap_or(0);
            let sum = left_strand + right_strand;
            if sum == 0 {
                [65535, 0]
            } else if sum == 1 {
                [0, 65535]
            } else {
                [0, 0]
            }
        })
        .collect();

    //TODO: fix ploidy_missingness, not correctly read from data
    let data_block = DataBlock {
        number_individuals,
        number_alleles: 2,
        minimum_ploidy: 2,
        maximum_ploidy: 2,
        ploidy_missingness,
        phased: false,
        bytes_probability: 16,
        probabilities,
    };

    //TODO: fix size_in_bytes, alleles if multiallelic ?
    let variant_data = VariantData {
        number_individuals: Some(number_individuals),
        variants_id: variant_id_fmt.to_string(),
        rsid: variant_id_fmt.to_string(),
        chr: chr.to_string(),
        pos: pos.parse().unwrap(),
        number_alleles: 2,
        alleles: vec![a1.to_string(), a2.to_string()],
        file_start_position: 0,
        size_in_bytes: 0,
        data_block,
    };
    Ok(variant_data)
}

fn parse_genotype_field(input: &str) -> IResult<&str, Vec<&str>> {
    let geno_start = "GT:AD:MD:DP:GQ:PL";
    // parse line until genotype starts
    let before_genotype_parser = preceded(preceded(take_until(geno_start), tag(geno_start)), tab);
    // parse genotype from list of values
    let parse_geno = terminated(take(3u8), take_while1(|c| c != '\t'));
    // parse whole line
    preceded(before_genotype_parser, separated_list0(tab, parse_geno))(input)
}

fn format_id_with_alleles(id: &str, a1: &str, a2: &str) -> Result<String> {
    let mut id_start: Vec<&str> = id.split(':').take(2).collect();
    id_start.push(a1);
    id_start.push(a2);
    Ok(id_start.join(":"))
}
