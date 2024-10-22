use bgen_reader::bgen::bgen_stream::write_samples;
use bgen_reader::bgen::header::{Header, HeaderFlags};
use bgen_reader::bgen::variant_data::{DataBlock, VariantData};
use clap::Parser;
use color_eyre::Report;
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

#[derive(Debug)]
enum VcfError {
    Io(std::io::Error),
    Nom(Report),
    Bgen(Report),
}

impl From<std::io::Error> for VcfError {
    fn from(error: std::io::Error) -> Self {
        VcfError::Io(error)
    }
}

impl From<Report> for VcfError {
    fn from(error: Report) -> Self {
        VcfError::Bgen(error)
    }
}

impl From<nom::Err<nom::error::Error<&str>>> for VcfError {
    fn from(error: nom::Err<nom::error::Error<&str>>) -> Self {
        VcfError::Nom(Report::msg(format!("Nom Error: {:?}", error)))
    }
}

#[derive(Parser, Debug)]
struct Args {
    /// Path to the input vcf file
    #[arg(short, long)]
    input: String,

    /// Path to the output bgen file
    #[arg(short, long)]
    output: String,
}

fn main() -> Result<(), VcfError> {
    let args = Args::parse();
    // let read_file = "/pasteur/zeus/projets/p02/GGS_WKD/HOME_LEO/vcf_to_bgen/data/small_2.vcf.gz";
    let mut reader = BufReader::new(MultiGzDecoder::new(File::open(&args.input)?));
    // First pass to get the number of variants
    let mut number_geno_line = 0;
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
            // If variant is multiallelic, we should add more than 1
            variant_num += alt_allele_count(&line)?;
            number_geno_line += 1;
        }
        line.clear();
    }
    bar.finish();
    println!("Done");
    drop(reader);

    // read file again
    let mut reader = BufReader::new(MultiGzDecoder::new(File::open(&args.input)?));

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
    let samples_str = samples_result?.1;
    let samples: Vec<String> = samples_str.into_iter().map(|s| s.to_string()).collect();
    let number_individuals = samples.len() as u32;
    line.clear();

    // write bgen file
    let mut bgen_writer = BufWriter::new(File::create(&args.output)?);
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
    let bar = ProgressBar::new(number_geno_line as u64);

    for _geno_line in 0..number_geno_line {
        reader.read_line(&mut line)?;
        let mut variant_data = parse_genotype_line(&line, number_individuals)?;
        let alt_variants: Vec<_> = variant_data.alleles[1]
            .split(',')
            .map(|s| s.to_string())
            .collect();
        if alt_variants.len() > 1 {
            // split multiallelic into biallelic
            for (alt_i, alt) in alt_variants.into_iter().enumerate() {
                let mut variant_data_clone = variant_data.clone();
                variant_data_clone.alleles[1] = alt.to_string();
                let variant_id_fmt = format_id_with_alleles(
                    &(variant_data.chr.to_string() + "_" + &variant_data.pos.to_string()),
                    &variant_data.alleles[0],
                    &alt,
                );
                variant_data_clone.variants_id = variant_id_fmt.clone();
                variant_data_clone.rsid = variant_id_fmt;
                // normalize probabilities
                variant_data_clone.data_block.probabilities = variant_data
                    .data_block
                    .probabilities
                    .chunks(2)
                    .flat_map(|g| genos_to_proba(g, alt_i as u32 + 1))
                    .collect();
                variant_data_clone.write_self(&mut bgen_writer, 2)?;
            }
        } else {
            // normalize probabilities
            variant_data.data_block.probabilities = variant_data
                .data_block
                .probabilities
                .chunks(2)
                .flat_map(|g| genos_to_proba(g, 1))
                .collect();
            variant_data.write_self(&mut bgen_writer, 2)?;
        }
        bar.inc(1);
        line.clear();
    }
    bar.finish();
    Ok(())
}

fn genos_to_proba(genos: &[u32], alt_number: u32) -> Vec<u32> {
    let left_strand = if genos[0] == alt_number { 1 } else { 0 };
    let right_strand = if genos[1] == alt_number { 1 } else { 0 };
    let sum = left_strand + right_strand;
    let result = if sum == 0 {
        [65535, 0]
    } else if sum == 1 {
        [0, 65535]
    } else {
        [0, 0]
    };
    result.to_vec()
}

fn parse_samples(input: &str) -> IResult<&str, Vec<&str>> {
    preceded(
        preceded(tag("#"), many0(preceded(alpha0, tab))),
        separated_list0(tab, alphanumeric0),
    )(input)
}

fn parse_one_field(input: &str) -> Result<(&str, &str), VcfError> {
    Ok(terminated(is_not("\t"), char('\t'))(input)?)
}

fn alt_allele_count(input: &str) -> Result<u32, VcfError> {
    let (remaining_input, _) = parse_one_field(input)?;
    let (remaining_input, _) = parse_one_field(remaining_input)?;
    let (remaining_input, _) = parse_one_field(remaining_input)?;
    let (remaining_input, _) = parse_one_field(remaining_input)?;
    let (_remaining_input, alt_alleles) = parse_one_field(remaining_input)?;
    Ok(alt_alleles.chars().filter(|&c| c == ',').count() as u32 + 1)
}

fn parse_genotype_line(input: &str, number_individuals: u32) -> Result<VariantData, VcfError> {
    let (remaining_input, chr) = parse_one_field(input)?;
    let (remaining_input, pos) = parse_one_field(remaining_input)?;
    let (remaining_input, variant_id) = parse_one_field(remaining_input)?;
    let (remaining_input, a1) = parse_one_field(remaining_input)?;
    let (remaining_input, a2) = parse_one_field(remaining_input)?;
    let genos_string = parse_genotype_field(remaining_input)?.1;
    let variant_id_fmt = format_id_with_alleles(variant_id, a1, a2);
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
            [left_strand, right_strand]
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

fn format_id_with_alleles(id: &str, a1: &str, a2: &str) -> String {
    let mut id_start: Vec<&str> = id.split(':').take(2).collect();
    id_start.push(a1);
    id_start.push(a2);
    id_start.join(":")
}
