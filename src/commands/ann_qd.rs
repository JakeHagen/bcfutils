use rust_htslib::bcf::{Read, Header};

pub fn ann_qd(input: Option<&str>, output: Option<&str>) {
    let mut bcf = bcfutils::get_rdr(input);
    
    let hdrv = bcf.header();
    let mut hdr = Header::from_template(&hdrv);
    hdr.push_record(r#"##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">"#.as_bytes());
    
    let mut obcf = bcfutils::get_wrtr(output, &hdr);

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
