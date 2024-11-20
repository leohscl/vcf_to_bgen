extern crate vcf_to_bgen;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader};
use vcf_to_bgen::{parse_genotype_line, read_vcf_header, split_multiallelic};

#[test]
fn read_samples() {
    let input = "data/only_header_10_samples.vcf.gz";
    // reads header
    let mut reader = BufReader::new(MultiGzDecoder::new(File::open(input).unwrap()));
    let samples = read_vcf_header(&mut reader).unwrap();
    assert_eq!(
        samples,
        [
            "HG00096", "HG00097", "HG00099", "HG00100", "HG00101", "HG00102", "HG00103", "HG00104",
            "HG00105", "HG001"
        ]
        .to_vec()
    );
}

#[test]
fn read_one_line() {
    let input = "data/100_vars_chr22_HG.vcf.gz";
    // reads header
    let mut reader = BufReader::new(MultiGzDecoder::new(File::open(input).unwrap()));
    let _samples = read_vcf_header(&mut reader).unwrap();
    // read first line
    let mut line = String::new();
    reader.read_line(&mut line).unwrap();
    let num_bits = 16;
    let number_individuals = 2504;
    let variant_data = parse_genotype_line(&line, number_individuals, num_bits).unwrap();
    let vec_variant_data = split_multiallelic(variant_data).unwrap();
    assert_eq!(
        vec_variant_data[0].data_block.probabilities[0..10],
        [65535, 0, 65535, 0, 65535, 0, 65535, 0, 65535, 0].to_vec()
    );
}

#[test]
fn read_one_line_2_field_format() {
    let input = "data/1_var_10_ind.vcf.gz";
    // reads header
    let mut reader = BufReader::new(MultiGzDecoder::new(File::open(input).unwrap()));
    let _samples = read_vcf_header(&mut reader).unwrap();
    // read first line
    let mut line = String::new();
    reader.read_line(&mut line).unwrap();
    let num_bits = 8;
    let number_individuals = 10;
    let variant_data = parse_genotype_line(&line, number_individuals, num_bits).unwrap();
    let vec_variant_data = split_multiallelic(variant_data).unwrap();
    assert_eq!(
        vec_variant_data[0].data_block.probabilities[0..10],
        [255, 0, 255, 0, 255, 0, 255, 0, 255, 0].to_vec()
    );
}

#[test]
fn read_one_line_complicated_format() {
    let input = "data/complicated_format.vcf.gz";
    // reads header
    let mut reader = BufReader::new(MultiGzDecoder::new(File::open(input).unwrap()));
    let _samples = read_vcf_header(&mut reader).unwrap();
    // read first line
    let mut line = String::new();
    reader.read_line(&mut line).unwrap();
    let num_bits = 8;
    let number_individuals = 10;
    let variant_data = parse_genotype_line(&line, number_individuals, num_bits).unwrap();
    let vec_variant_data = split_multiallelic(variant_data).unwrap();
    assert_eq!(
        vec_variant_data[0].data_block.probabilities[0..10],
        [255, 0, 255, 0, 255, 0, 255, 0, 255, 0].to_vec()
    );
}

#[test]
fn read_one_line_missing_values() {
    let input = "data/1_var_10_ind_with_missing.vcf.gz";
    // reads header
    let mut reader = BufReader::new(MultiGzDecoder::new(File::open(input).unwrap()));
    let _samples = read_vcf_header(&mut reader).unwrap();
    // read first line
    let mut line = String::new();
    reader.read_line(&mut line).unwrap();
    let num_bits = 8;
    let number_individuals = 10;
    let variant_data = parse_genotype_line(&line, number_individuals, num_bits).unwrap();
    let vec_variant_data = split_multiallelic(variant_data).unwrap();
    // probabilities are not impacted by missing values
    assert_eq!(
        vec_variant_data[0].data_block.probabilities[0..10],
        [255, 0, 255, 0, 255, 0, 255, 0, 255, 0].to_vec()
    );
    assert_eq!(
        vec_variant_data[0].data_block.ploidy_missingness[0..10],
        [130, 2, 130, 130, 2, 2, 2, 2, 2, 2].to_vec()
    );
}

#[test]
fn read_one_line_multiallelic() {
    let input = "data/multiallelic_1_var.vcf.gz";
    // reads header
    let mut reader = BufReader::new(MultiGzDecoder::new(File::open(input).unwrap()));
    let _samples = read_vcf_header(&mut reader).unwrap();
    // read first line
    let mut line = String::new();
    reader.read_line(&mut line).unwrap();
    let num_bits = 8;
    let number_individuals = 10;
    let variant_data = parse_genotype_line(&line, number_individuals, num_bits).unwrap();
    let vec_variant_data = split_multiallelic(variant_data).unwrap();
    assert_eq!(
        vec_variant_data[0].data_block.probabilities[0..10],
        vec![255, 0, 255, 0, 0, 255, 255, 0, 255, 0]
    );
    assert_eq!(
        vec_variant_data[0].data_block.ploidy_missingness[0..10],
        [2, 130, 2, 2, 130, 2, 2, 2, 2, 2].to_vec()
    );
    assert_eq!(
        vec_variant_data[1].data_block.probabilities[0..10],
        vec![255, 0, 0, 255, 255, 0, 255, 0, 0, 255]
    );
    assert_eq!(
        vec_variant_data[1].data_block.ploidy_missingness[0..10],
        [2, 2, 130, 2, 2, 2, 2, 130, 2, 2].to_vec()
    );
}
