FROM rust:latest as builder

WORKDIR /usr/src/vcf_to_bgen
COPY . .
COPY ./bgen_reader/ ../bgen_reader
RUN cargo install --path .

FROM debian:latest
RUN apt-get update && apt-get -y install  sqlite3 && rm -rf /var/lib/apt/lists/*
COPY --from=builder /usr/local/cargo/bin/vcf_to_bgen /usr/local/bin/vcf_to_bgen
CMD ["vcf_to_bgen"]
