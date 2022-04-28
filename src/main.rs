extern crate rust_htslib;
extern crate serde;
extern crate bio;
extern crate linear_map;
use bio::io::gff;
use crate::rust_htslib::bcf::{Reader, Read, Writer, Format, Header};
use rust_htslib::bcf::header::HeaderRecord;
use rust_htslib::bcf::record::{Buffer};
use linear_map::LinearMap;
use serde::Deserialize;
use std::collections::HashMap;
use std::error::Error;
use std::str;
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
    },
    MCSQ {
        input: Option<String>,
        #[clap(long, short)]
        output: Option<String>,
        #[clap(long, short)]
        gff: Option<String>,
    },
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

//fn extract_transcript() -> {
//
//}
//
fn build_trx_map(gff_fp: Option<&str>) -> HashMap<String, String> {
    let mut gff_rdr = gff::Reader::from_file(gff_fp.unwrap(), gff::GffType::GFF3).expect("error reading gff");
    let mut trx_map: HashMap<String, String> = HashMap::new();
    for record in gff_rdr.records() {
        let mut s = String::new();
        s.push('|');
        let rec = record.ok().expect("Error reading record.");
        if rec.feature_type() == "transcript" {
            let tags = rec.attributes().get_vec("tag");
            match tags {
                None => continue,
                Some(v) => {
                    if v.contains(&"Ensembl_canonical".to_string()) {
                        s.push_str(&"YES");
                    }
                }
            }
            trx_map.insert(rec.attributes().get("ID").expect("transcript did not have ID for some reason").split('.').next().unwrap().to_string().clone(), s);
        }
    }
    return trx_map;
}

fn get_bcsq_hdr_map(hdr_recs: Vec<HeaderRecord>) -> Option<LinearMap<String, String>> { //rust_htslib::bcf::HeaderRecord> {
    for hrec in hdr_recs.iter() {
        match hrec {
            HeaderRecord::Info{key,values} => {
                if values.get("ID").unwrap() == "BCSQ" {
                    return Some(values.clone())
                }
            }
            _ => continue,
        }
    }
    return None
}

fn get_num_bcsq_keys(desc: &String) -> usize {
    return desc.split('|').count()
}


fn mcsq(input: Option<&str>, output: Option<&str>, gff_fp: Option<&str>) {
    let mut bcf = get_rdr(input);
    let hdrv = bcf.header();

    let bcsq_map = get_bcsq_hdr_map(hdrv.header_records()).expect("was not able to get BCSQ header info line");
    let num_keys = get_num_bcsq_keys(bcsq_map.get("Description").expect("bcsq map doesnt have \"Description\", it really should though, something funky is happening"));
    
    let mut hdr = Header::from_template(&hdrv);
    hdr.remove_info(b"BCSQ");
    hdr.push_record(r#"##INFO=<ID=BCSQ,Number=.,Type=String,Description="Local consequence annotation from BCFtools/csq, see http://samtools.github.io/bcftools/howtos/csq-calling.html for details. Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change|CANONICAL">"#.as_bytes());
    let mut obcf = get_wrtr(output, &hdr);

    let trx_map = build_trx_map(gff_fp);
    let mut b = Buffer::new();
    for record_result in bcf.records() {
        let mut record = record_result.expect("fail to read record");
        let bcsqs = record.info_shared_buffer(b"BCSQ", &mut b).string().unwrap().unwrap();
        let mut mbcsqs = String::new();
        for (i, bcsq_b) in bcsqs.iter().enumerate() {
            let num_c_bcsq = get_num_bcsq_keys(bcsq_b.to_string()); 
            if i > 0 {
                mbcsqs.push(',');
            };
            let bcsq = str::from_utf8(bcsq_b).unwrap();
            mbcsqs.push_str(bcsq);
            let trn = bcsq.split('|').collect::<Vec<&str>>()[2];
            match trx_map.get(trn) {
                Some(t) => mbcsqs.push_str(t),
                None => mbcsqs.push('|'),
            }
        };
        record.push_info_string(b"BCSQ", &[mbcsqs.as_bytes()]).expect("failed to set QD info field");
        obcf.write(&record).expect("failed to write record");

        //let bcsq = match str::from_utf8(record.info(b"BCSQ").string().unwrap()) {
        //    Ok(v) => v,
        //    Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        //};
        //println!("{}", bcsq);
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
        Commands::MCSQ { input, output, gff } => {
            mcsq(input.as_deref(), output.as_deref(), gff.as_deref())
        }
    }
}
