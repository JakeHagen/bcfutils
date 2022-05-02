use std::collections::HashMap;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::str;
use crate::rust_htslib::bcf::{Read, Header};
use rust_htslib::bcf::header::HeaderRecord;
use rust_htslib::bcf::record::Buffer;
use linear_map::LinearMap;

fn extract_trns_id(s: &str) -> &str {
    let sb = s.find("ID=").expect("transcript line doesnt have ID= field");
    let peb = sb + s[sb..].find(";").expect("could not find next ';', line is malformed");
    let eb = sb + s[sb..peb].find(".").expect("transcript doesnt have '.' notation");
    return &s[sb+3..eb];
}

fn get_field<'a>(s: &'a str, field: &'a str) -> Option<&'a str> {
    let sb = match s.find(field) {
        Some(b) => b+field.len()+1, // +1 for the equal sign
        None => return None,
    };
    let eb = match s[sb..].find(";") {
        Some(b) => sb + b,
        None => s.len(),
    };
    return Some(&s[sb..eb]);
}

fn has_tag(tags: &str, tag: &str) -> bool {
    match tags.find(tag) {
        Some(_) => return true,
        None => return false,
    };
}

fn extract_appris(tags: &str) -> Option<&str> {
    let sb = match tags.find("appris_") {
        Some(b) => b+7, // 7 == length of appris_
        None => return None,
    };
    let eb = match tags[sb..].find(",") {
        Some(b) => sb + b,
        None => tags.len(),
    };
    return Some(&tags[sb..eb]);
}

fn extract_uncertain_start_end(tags: &str) -> String {
    let mut found_tags = vec![];
    let possible_tags = vec!["cds_start_NF", "cds_end_NF", "mRNA_start_NF", "mRNA_end_NF"];
    for tag in possible_tags {
        if has_tag(tags, tag) {
            found_tags.push(tag);
        }
    }
    return found_tags.join("&");
}

fn add_trns_to_map(line: String, map: &mut HashMap<String, String>) {
    for (i, col) in line.split("\t").enumerate() {
        if i == 2 {
            if col != "transcript" { return; }
        }
        if i == 8 {
            let mut s = String::new();
            match get_field(col, "tag") {
                Some(tags) => {
                    // add CANONICAL
                    if has_tag(tags, "Ensembl_canonical") {
                        s.push_str("|YES")
                    } else {
                        s.push('|')
                    }
                    // add appris
                    match extract_appris(tags) {
                        Some(appris) => {
                            s.push('|');
                            s.push_str(appris);
                        }
                        None => s.push('|'),
                    }
                    // add CCDS
                    if has_tag(tags, "CCDS") {
                        s.push_str("|CCDS")
                    } else {
                        s.push('|')
                    }
                    // add uncertain start/end
                    s.push('|');
                    s.push_str(&extract_uncertain_start_end(tags));
                }
                None => s.push_str("||||"),
            };
            // add TSL
            match get_field(col, "transcript_support_level") {
                Some(tsl) => {
                    s.push('|');
                    s.push_str(tsl);
                },
                None => s.push('|'),
            };

            let trns_id = extract_trns_id(col);
            map.insert(trns_id.to_string(), s);
        }
    }
}

fn build_trx_map(gff_fp: Option<&str>) -> HashMap<String, String> {
    let gff = File::open(gff_fp.unwrap());
    let rdr = BufReader::new(gff.expect("couldnt make buffer for gff file"));

    let mut trns_map: HashMap<String, String> = HashMap::new();

    for line in rdr.lines() {
        if let Ok(l) = line {
            if &l[0..1] == "#" {
                continue;
            }
            add_trns_to_map(l, &mut trns_map);
        }
    }
    return trns_map;
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

fn get_num_bcsq_keys(desc: &str) -> usize {
    return desc.split('|').count()
}


pub fn mcsq(input: Option<&str>, output: Option<&str>, gff_fp: Option<&str>) {
    let mut bcf = bcfutils::get_rdr(input);
    let hdrv = bcf.header();

    let bcsq_map = get_bcsq_hdr_map(hdrv.header_records()).expect("was not able to get BCSQ header info line");
    let num_keys = get_num_bcsq_keys(bcsq_map.get("Description").expect("bcsq map doesnt have \"Description\", it really should though, something funky is happening"));

    let mut hdr = Header::from_template(&hdrv);
    hdr.remove_info(b"BCSQ");
    hdr.push_record(r#"##INFO=<ID=BCSQ,Number=.,Type=String,Description="Local consequence annotation from BCFtools/csq, see http://samtools.github.io/bcftools/howtos/csq-calling.  html for details. Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change|CANONICAL|appris|ccds|unknown_start_end|TSL">"#.as_bytes());
    let mut obcf = bcfutils::get_wrtr(output, &hdr);

    let trx_map = build_trx_map(gff_fp);
    let mut b = Buffer::new();
    for record_result in bcf.records() {
        let mut record = record_result.expect("fail to read record");
        let bcsqs = record.info_shared_buffer(b"BCSQ", &mut b).string().unwrap().unwrap();
        let mut mbcsqs = String::new();
        for (i, bcsq_b) in bcsqs.iter().enumerate() {
            if i > 0 {
                mbcsqs.push(',');
            };
            let bcsq = str::from_utf8(bcsq_b).unwrap();
            let num_c_bcsq = get_num_bcsq_keys(bcsq);

            mbcsqs.push_str(bcsq);

            let num_to_fill = num_keys - num_c_bcsq;
            for _ in 0..num_to_fill {
                mbcsqs.push('|');
            }
            let trn = bcsq.split('|').collect::<Vec<&str>>()[2];
            match trx_map.get(trn) {
                Some(t) => mbcsqs.push_str(t),
                None => mbcsqs.push_str("|||||"),
            }
        };
        record.push_info_string(b"BCSQ", &[mbcsqs.as_bytes()]).expect("failed to set QD info field");
        obcf.write(&record).expect("failed to write record");
    }
}


