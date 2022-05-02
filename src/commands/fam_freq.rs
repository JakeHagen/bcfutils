use csv;
use std::error::Error;
use serde::Deserialize;

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

pub fn fam_freq(pedigree: Option<&str>) {
    let r = parse_ped(pedigree);
}
