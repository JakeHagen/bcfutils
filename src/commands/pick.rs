use std::str;
use crate::rust_htslib::bcf::{Read, Header};
use rust_htslib::bcf::record::Buffer;
use std::collections::HashMap;

//fn spllit_csqs<'a>(csqs: String) -> Vec<&'a str> {
//    return csqs.split(",").collect::<Vec<&str>>()
//}

fn extract_field<'a>(csq: &'a str, index: usize) -> &'a str {
    return csq.split("|").collect::<Vec<&str>>()[index]
}

fn get_canon_rank(csq: &str) -> u8 {
    let canon = extract_field(csq, 7);
    return match canon {
        "YES" => 1,
        _ => 99,
    }
}

fn get_appris_rank(appris: &str) -> u8 {
    return match appris {
        "principle_1" => 1,
        "principle_2" => 2,
        "principle_3" => 3,
        "principle_4" => 4,
        "principle_5" => 5,
        "alternative_2" => 6,
        "alternative_1" => 7,
        _ => 99,
    }
}

fn get_tsl_rank(csq: &str) -> u8 {
    let tsl = extract_field(csq, 11);
    return match tsl {
        "1" => 1,
        "2" => 2, 
        "3" => 3,
        "4" => 4,
        "5" => 5,
        _ => 99
    }
}


fn csq_filter<'a>(csqs: Vec<&'a str>, get_rank: &'a fn(&str) -> u8) -> Vec<&'a str> {
    let mut f_csqs: Vec<&str> = vec![];
    let mut current_min = 99;
    for csq in csqs.iter() {
        let rank = get_rank(csq);
        if rank <= current_min {
            f_csqs.push(csq);
            current_min = rank;
        }
    }
    if f_csqs.len() > 0 {
        return f_csqs;
    }
    return csqs;
}

//, &get_appris_rank, &get_tsl_rank
fn pick_one(csqs: Vec<&str>) -> &str {
    let criteria_funcs: Vec<&dyn Fn(&str) -> u8> = vec![&get_canon_rank, &get_appris_rank, &get_tsl_rank];
    for cf in criteria_funcs {
       if csqs.len() == 1 {
            return csqs[0];
        }
        let csqs = csq_filter(csqs, cf);
    }
    return csqs[0];
}

//Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change|CANONICAL|appris|ccds|unknown_start_end|TSL
pub fn pick(input: Option<&str>, output: Option<&str>) {
    let mut bcf = bcfutils::get_rdr(input);
    let hdrv = bcf.header();

    let mut hdr = Header::from_template(&hdrv);
    hdr.push_record(r#"##INFO=<ID=pick_Consequence,Number=1,Type=String,Description="picked one consequence from BCFtools/csq + mcsq">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=pick_gene,Number=1,Type=String,Description="picked one consequence from BCFtools/csq + mcsq">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=pick_transcript,Number=1,Type=String,Description="picked one consequence from BCFtools/csq + mcsq">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=pick_biotype,Number=1,Type=String,Description="picked one consequence from BCFtools/csq + mcsq">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=pick_strand,Number=1,Type=String,Description="picked one consequence from BCFtools/csq + mcsq">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=pick_amino_acid_change,Number=1,Type=String,Description="picked one consequence from BCFtools/csq + mcsq">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=pick_dna_change,Number=1,Type=String,Description="picked one consequence from BCFtools/csq + mcsq">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=pick_CANONICAL,Number=1,Type=String,Description="picked one consequence from BCFtools/csq + mcsq">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=pick_appris,Number=1,Type=String,Description="picked one consequence from BCFtools/csq + mcsq">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=pick_ccds,Number=1,Type=String,Description="picked one consequence from BCFtools/csq + mcsq">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=pick_unknown_start_end,Number=1,Type=String,Description="picked one consequence from BCFtools/csq + mcsq">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=pick_TSL,Number=1,Type=String,Description="picked one consequence from BCFtools/csq + mcsq">"#.as_bytes());
    let mut obcf = bcfutils::get_wrtr(output, &hdr);

    let mut b = Buffer::new();
    for record_result in bcf.records() {
        let mut record = record_result.expect("fail to read record");
        let bcsqs = match record.info_shared_buffer(b"BCSQ", &mut b).string().unwrap() {
            Some(b) => b,
            None => {
                obcf.write(&record).expect("failed to write record");
                continue;
            }
        };

        let mut bcsqs_v = vec![];
        for bcsq_b in bcsqs.iter() {
            let bcsq = str::from_utf8(bcsq_b).unwrap();
            bcsqs_v.push(bcsq);
            //if canon(bcsq) {
            //    bcsqs_v.push(bcsq);
            //}
        };


        let picked = pick_one(bcsqs_v);
        record.push_info_string(b"pick_Consequence", &[picked.as_bytes()]).expect("failed to set pick_Consequence field");
        obcf.write(&record).expect("failed to write record");
    }
}
