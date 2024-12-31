use bgen_reader::bgen::bgen_stream::write_samples;
use bgen_reader::bgen::header::{Header, HeaderFlags};
use bgen_reader::bgen::variant_data::{DataBlock, VariantData};
use color_eyre::Report;
use flate2::read::MultiGzDecoder;
use indicatif::ProgressBar;
use nom::branch::alt;
use nom::bytes::complete::{is_not, tag, take, take_while1};
use nom::character::complete::{alpha0, alphanumeric0, char, tab};
use nom::combinator::success;
use nom::multi::{count, many0, separated_list0};
use nom::sequence::{delimited, preceded, terminated};
use nom::{IResult, InputIter};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::time::Duration;

#[derive(Debug)]
pub enum VcfError {
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

// Wrapper type for variant data, with added genotype represented as vcf string
pub struct VariantDataToParse<'a> {
    variant_data: VariantData,
    geno_string_vcf: Vec<&'a str>,
}

pub fn count_variants(input: &str) -> Result<(u32, u32), VcfError> {
    let mut reader = BufReader::new(MultiGzDecoder::new(File::open(input)?));
    let mut number_geno_line = 0;
    let mut variant_num = 0;
    let mut line = String::new();
    println!("Counting variants...  ");
    let bar = ProgressBar::new_spinner();
    bar.enable_steady_tick(Duration::from_millis(100));
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
    Ok((variant_num, number_geno_line))
}

pub fn read_vcf_header(reader: &mut impl BufRead) -> Result<Vec<String>, VcfError> {
    let mut line = String::new();
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
    Ok(samples_str.into_iter().map(|s| s.to_string()).collect())
}

pub fn write_bgen_header(
    bgen_writer: &mut BufWriter<std::fs::File>,
    samples: &[String],
    number_individuals: u32,
    variant_num: u32,
) -> Result<(), VcfError> {
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
    header.write_header(bgen_writer)?;
    //
    // write samples
    Ok(write_samples(samples, bgen_writer, len_sample_block)?)
}

pub fn parse_geno_line(
    vec_probas: &mut [u32],
    vec_ploidy_m: &mut [u8],
    geno_line: &[&str],
    alt_allele_num: usize,
    num_bits: u8,
) {
    geno_line.iter().enumerate().for_each(|(geno_i, geno_s)| {
        let mut geno_iter = geno_s
            .iter_elements()
            .filter_map(|c| c.to_digit(10))
            .filter(|&d| d == 0 || d == alt_allele_num as u32)
            .map(|d| if d == 0 { 0 } else { 1 });
        let count_valid = geno_iter.clone().count();
        // if there is less than 2 values, there is missingness
        let ploidy_m = if count_valid < 2 { (1u8 << 7) + 2 } else { 2u8 };
        let left_strand = geno_iter.next().unwrap_or(0);
        let right_strand = geno_iter.next().unwrap_or(0);
        let genos = [left_strand, right_strand];
        // convert geno to bgen probabilities
        let probas = genos_to_proba(&genos, num_bits);
        vec_probas[geno_i * 2] = probas[0];
        vec_probas[geno_i * 2 + 1] = probas[1];
        vec_ploidy_m[geno_i] = ploidy_m;
    });
}

pub fn parse_vcf_geno(
    variant_data_to_parse: &VariantDataToParse<'_>,
    alt_allele: String,
    alt_allele_num: usize,
    num_bits: u8,
    number_individuals: u32,
) -> VariantData {
    let number_individuals = number_individuals as usize;
    // use variant data as pattern
    let mut variant_data_clone = variant_data_to_parse.variant_data.clone();

    // fill description fields
    let variant_id_fmt = format_id_with_alleles(
        &(variant_data_clone.chr.to_string() + ":" + &variant_data_clone.pos.to_string()),
        &variant_data_clone.alleles[0],
        &alt_allele,
    );
    variant_data_clone.variants_id = variant_id_fmt.clone();
    variant_data_clone.alleles[1] = alt_allele;
    variant_data_clone.rsid = variant_id_fmt;

    let mut ploidy_missingness = vec![0; number_individuals];
    let mut probabilities = vec![0; number_individuals * 2];

    // convert string to missingness and probas
    parse_geno_line(
        &mut probabilities,
        &mut ploidy_missingness,
        &variant_data_to_parse.geno_string_vcf,
        alt_allele_num,
        num_bits,
    );
    variant_data_clone.data_block.ploidy_missingness = ploidy_missingness;
    variant_data_clone.data_block.probabilities = probabilities;
    variant_data_clone
}

pub fn split_multiallelic(
    variant_data_to_parse: VariantDataToParse<'_>,
    number_individuals: u32
) -> Result<Vec<VariantData>, VcfError> {
    let variant_data = &variant_data_to_parse.variant_data;

    let alt_variants: Vec<_> = variant_data.alleles[1]
        .split(',')
        .map(|s| s.to_string())
        .collect();
    let num_bits = variant_data.data_block.bits_storage;
    // split multiallelic into biallelic
    let vec_variant_data = alt_variants
        .into_iter()
        .enumerate()
        .map(|(alt_i, alt)| parse_vcf_geno(&variant_data_to_parse, alt, alt_i + 1, num_bits, number_individuals))
        .collect::<Vec<VariantData>>();
    Ok(vec_variant_data)
}

pub fn convert_variant_blocks(
    reader: &mut impl BufRead,
    bgen_writer: &mut BufWriter<std::fs::File>,
    number_geno_line: u32,
    number_individuals: u32,
    num_bits: u8,
) -> Result<(), VcfError> {
    let mut line = String::new();

    let bar = ProgressBar::new(number_geno_line as u64);

    for _geno_line in 0..number_geno_line {
        reader.read_line(&mut line)?;
        let variant_data = parse_genotype_line(&line, number_individuals, num_bits)?;
        let vec_variant_data = split_multiallelic(variant_data, number_individuals)?;
        for var_data in vec_variant_data {
            var_data.write_self(bgen_writer, 2)?;
        }
        bar.inc(1);
        line.clear();
    }
    bar.finish();
    Ok(())
}

pub fn convert_to_bgen(
    input: &str,
    output: &str,
    variant_num: u32,
    number_geno_line: u32,
    num_bits: u8,
) -> Result<(), VcfError> {
    // reads vcf
    let mut reader = BufReader::new(MultiGzDecoder::new(File::open(input)?));
    // writes bgen
    let mut bgen_writer = BufWriter::new(File::create(output)?);

    // get samples from header
    let samples = read_vcf_header(&mut reader)?;
    let number_individuals = samples.len() as u32;

    // write header and samples
    write_bgen_header(&mut bgen_writer, &samples, number_individuals, variant_num)?;

    // write variant blocks
    println!("Converting variants to bgen format");
    convert_variant_blocks(
        &mut reader,
        &mut bgen_writer,
        number_geno_line,
        number_individuals,
        num_bits,
    )
}

fn genos_to_proba(genos: &[u32], num_bits: u8) -> Vec<u32> {
    let sum = genos[0] + genos[1];
    let proba_1 = (1 << num_bits) - 1;
    let result = if sum == 0 {
        [proba_1, 0]
    } else if sum == 1 {
        [0, proba_1]
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

pub fn parse_genotype_line(
    input: &str,
    number_individuals: u32,
    num_bits: u8,
) -> Result<VariantDataToParse<'_>, VcfError> {
    let (remaining_input, chr) = parse_one_field(input)?;
    let (remaining_input, pos) = parse_one_field(remaining_input)?;
    let (remaining_input, variant_id) = parse_one_field(remaining_input)?;
    let (remaining_input, a1) = parse_one_field(remaining_input)?;
    let (remaining_input, a2) = parse_one_field(remaining_input)?;
    let genos_string = parse_genotype_field(remaining_input)?.1;
    let variant_id_fmt = format_id_with_alleles(variant_id, a1, a2);
    let data_block = DataBlock {
        number_individuals,
        number_alleles: 2,
        minimum_ploidy: 2,
        maximum_ploidy: 2,
        ploidy_missingness: vec![],
        phased: false,
        bits_storage: num_bits,
        probabilities: vec![],
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
    let variant_data_to_parse = VariantDataToParse {
        variant_data,
        geno_string_vcf: genos_string,
    };
    Ok(variant_data_to_parse)
}

fn parser_elt_tab(input: &str) -> IResult<&str, &str> {
    let until_tab = take_while1(|c| c != '\t');
    terminated(until_tab, tab)(input)
}

fn parser_elt_colon(input: &str) -> IResult<&str, &str> {
    terminated(is_not(":"), tag(":"))(input)
}

fn parse_genotype_field(input: &str) -> IResult<&str, Vec<&str>> {
    //// V1
    //let geno_start = "GT:AD:MD:DP:GQ:PL";
    //// parse line until genotype starts
    //let before_genotype_parser = preceded(preceded(take_until(geno_start), tag(geno_start)), tab);
    //// parse genotype from list of values
    //let parse_geno = terminated(take(3u8), take_while1(|c| c != '\t'));
    //// parse whole line
    //preceded(before_genotype_parser, separated_list0(tab, parse_geno))(input)

    //// V2
    //let geno_start = "GT:AD:MD:DP:GQ:PL";
    //// parse line until genotype starts
    //let before_genotype_parser = preceded(preceded(take_until(geno_start), tag(geno_start)), tab);
    //// parse genotype from list of values
    //let parse_geno = terminated(take(3u8), is_not("\t"));
    //// parse whole line
    //preceded(before_genotype_parser, separated_list0(tab, parse_geno))(input)


    // V3
    let until_tab = take_while1(|c| c != '\t');
    // Genotype starts at column 9, 5 lines are already read
    let mut before_genotype_parser = preceded(count(parser_elt_tab, 3), parser_elt_tab);
    // Gives Format field, and remaining line is left to parse
    let parse_line_start = before_genotype_parser(input).unwrap();
    // Format like GT:GP..
    let remaining_string = parse_line_start.0;
    let format = parse_line_start.1;
    let gt_position = format.split(':').position(|s| s == "GT").unwrap();

    // let parse_geno = delimited(count(parser_elt_colon, gt_position), take(3u8), is_not("\t"));
    let parse_geno = delimited(
        count(parser_elt_colon, gt_position),
        take(3u8),
        alt((until_tab, success("1"))),
    );
    separated_list0(tab, parse_geno)(remaining_string)
}

fn format_id_with_alleles(id: &str, a1: &str, a2: &str) -> String {
    let mut id_start: Vec<&str> = id.split(':').take(2).collect();
    id_start.push(a1);
    id_start.push(a2);
    id_start.join(":")
}
