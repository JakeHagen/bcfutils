extern crate rust_htslib;
extern crate serde;
use crate::rust_htslib::bcf::{Reader, Read, Writer, Format, Header};
use serde::Deserialize;
use std::env;
use std::io;
use std::process;
use std::error::Error;
use csv;
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    QD {
        input: Option<String>,
        #[clap(long, short)]
        output: Option<String>,
    },
    FamFreq {
        input: Option<String>,
        #[clap(long, short)]
        output: Option<String>,
        #[clap(long, short)]
        pedigree: Option<String>,
    }
}

fn get_rdr(input: Option<&str>) -> rust_htslib::bcf::Reader {
    match input {
        None => {
           return  Reader::from_stdin().expect("failed to create reader from stdin");
        }
        Some("-") => {
            return Reader::from_stdin().expect("failed to create reader from stdin");
        }
        Some(inner) => {
            return Reader::from_path(inner).expect("failed to create reader from file");
        }
    }
}

fn get_wrtr(input: Option<&str>, hdr: &rust_htslib::bcf::Header) -> rust_htslib::bcf::Writer {
    match input {
        None => {
            return Writer::from_stdout(hdr, false, Format::Bcf).expect("unable to create writer on stdout");
        }
        Some("-") => {
            return Writer::from_stdout(hdr, false, Format::Bcf).expect("unable to create writer on stdout");
        }
        Some(inner) => {
            return Writer::from_path(inner, hdr, false, Format::Bcf).expect("unable to create writer with filepath");
        }
    }
}

#[derive(Debug, Deserialize)]
struct Sample {
    fid: String,
    sid: String,
    pid: String,
    mid: String,
    sex: String,
    status: String,
}

fn parse_ped(pedigree: Option<&str>) -> Result<(), Box<dyn Error>> {
    let mut ped = match csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(pedigree.unwrap()) {
        Ok(ped) => ped,
        Err(error) => {
            panic!("failed to create reader: {}", error);
        }
    };
    
    for result in ped.deserialize() {
        let sample: Sample = result?;
    };
    Ok(())
}

fn fam_freq(pedigree: Option<&str>) {
    let r = parse_ped(pedigree);

}


fn ann_qd(input: Option<&str>, output: Option<&str>) {
    let mut bcf = get_rdr(input);
    let hdrv = bcf.header();
    let mut hdr = Header::from_template(&hdrv);
    hdr.push_record(r#"##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">"#.as_bytes());
    let mut obcf = get_wrtr(output, &hdr);

    for record_result in bcf.records() {
        let mut record = record_result.expect("fail to read record");
        obcf.translate(&mut record);
        let qual = record.qual();
        let mut dp = 0i32;
        let sample_count = usize::try_from(record.sample_count()).unwrap();
        let gts = record.genotypes().expect("Error reading genotypes");
        let dps = record.format(b"DP").integer().expect("Couldn't retrieve DP field");
        'sample_: for sidx in 0..sample_count {
            for gta in gts.get(sidx).iter() {
                match gta.index() {
                    None => continue 'sample_,
                    Some(0) => continue,
                    Some(_) => {
                        dp += dps[sidx][0];
                        continue 'sample_;
                    }
                }
            }
        }
        if dp > 0 {
            let qd = qual / dp as f32;
            record.push_info_float(b"QD", &[qd]).expect("failed to set QD info field");
        }
        obcf.write(&record).expect("failed to write record");
        
    }
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Commands::QD { input, output } => {
            ann_qd(input.as_deref(), output.as_deref())
        }
        Commands::FamFreq { input, output, pedigree } => {
            fam_freq(pedigree.as_deref())
        }
    }
}
