extern crate vcf_to_bgen;
use vcf_to_bgen::count_variants;

#[test]
fn count_100_variants() {
    let input = "data/100_vars_chr22_HG.vcf.gz";
    let (num_variant, num_geno_line) = count_variants(input).unwrap();
    assert_eq!(num_geno_line, 100);
    assert_eq!(num_variant, 100);
}

#[test]
fn count_variants_with_multiallelic() {
    let input = "data/10_vars_with_multiallelic.vcf.gz";
    let (num_variant, num_geno_line) = count_variants(input).unwrap();
    assert_eq!(num_geno_line, 10);
    assert_eq!(num_variant, 11);
}
