use csv;
use std::error::Error;
use std::collections::HashMap;
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

fn mk_sid_fid_map(pedigree: Option<&str>) -> Result<HashMap<String, String>, Box<dyn Error>> {
    let ped_path = match pedigree {
        None => panic!("pedigree must be specified with the \"-p\" option"),
        Some(ped) => ped
    };
    let mut ped = match csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(ped_path) {
        Ok(p) => p,
        Err(error) => {
            panic!("failed to create reader from {}: {}", ped_path, error);
        }
    };

    let mut map: HashMap<String, String> = HashMap::new();
    for result in ped.deserialize() {
        let sample: Sample = result?;
        map.insert(sample.sid, sample.fid);
    };
    return Ok(map)
}

pub fn fam_freq(pedigree: Option<&str>) {
    let map = mk_sid_fid_map(pedigree).expect("unable to make sample to family map");

}
