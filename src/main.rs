use clap::Parser;
use vcf_to_bgen::{convert_to_bgen, count_variants, VcfError};

#[derive(Parser, Debug)]
struct Args {
    /// Path to the input vcf file
    #[arg(short, long)]
    input: String,

    /// Path to the output bgen file
    #[arg(short, long)]
    output: String,

    /// Number of bits used for probability storage
    #[arg(long)]
    num_bits: Option<u8>,
}

fn main() -> Result<(), VcfError> {
    let args = Args::parse();
    // First pass to get the number of variants
    let (variant_num, number_geno_line) = count_variants(&args.input)?;
    // Convert to bgen, line by line
    convert_to_bgen(
        &args.input,
        &args.output,
        variant_num,
        number_geno_line,
        args.num_bits.unwrap_or(8),
    )?;
    Ok(())
}
