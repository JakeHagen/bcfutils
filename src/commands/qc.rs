use rust_htslib::bcf::Read;
use std::io;
use std::str;
use phf::phf_map;
use std::fs::File;

static TI: phf::Map<&[u8], phf::Map<&[u8], bool>> = phf_map! {
    b"A" => phf_map! {
        b"G" => true,
        b"T" => false,
        b"C" => false,
    },
    b"C" => phf_map! {
        b"T" => true,
        b"A" => false,
        b"G" => false,
    },
    b"G" => phf_map! {
        b"A" => true,
        b"C" => false,
        b"T" => false,
    },
    b"T" => phf_map! {
        b"C" => true,
        b"A" => false,
        b"G" => false,
    },
};

struct Metrics {
    n_vars: Vec<i32>,
    n_indels: Vec<i32>,
    n_snvs: Vec<i32>,
    n_hets: Vec<i32>,
    n_hom_alts: Vec<i32>,
    n_tis: Vec<i32>,
    n_tvs: Vec<i32>,
}

impl Metrics {
    fn new(n_samples: usize) -> Metrics {
        return Metrics {
            n_vars: vec![0; n_samples],
            n_indels: vec![0; n_samples],
            n_snvs: vec![0; n_samples],
            n_hets: vec![0; n_samples],
            n_hom_alts: vec![0; n_samples],
            n_tis: vec![0; n_samples],
            n_tvs: vec![0; n_samples],
        }
    }

    fn update(&mut self, sidx: usize, gt: &[i32], alleles: &Vec<&[u8]>) {
        self.n_vars[sidx] += 1;
        if is_hom(gt) {
            self.n_hom_alts[sidx] += 1;
            if is_snv(decode_genotype(gt[0]), &alleles) {
                self.n_snvs[sidx] += 1;
                if is_ti(alleles[0], alleles[decode_genotype(gt[1]) as usize]) {
                    self.n_tis[sidx] += 1
                } else {
                    self.n_tvs[sidx] += 1
                }
            } else {
                self.n_indels[sidx] += 1
            }
        } else {
            self.n_hets[sidx] += 1;
            for g in gt {
                if decode_genotype(*g) <= 0 {
                    continue;
                }
                if is_snv(decode_genotype(*g), &alleles) {
                    self.n_snvs[sidx] += 1;
                    if is_ti(alleles[0], alleles[decode_genotype(*g) as usize]) {
                        self.n_tis[sidx] += 1
                    } else {
                        self.n_tvs[sidx] += 1
                    }
                } else {
                    self.n_indels[sidx] += 1
                }
            }
        }

    }
}

fn decode_genotype(g: i32) -> i32 {
    return (g >> 1) - 1
}

fn is_ti(r: &[u8], a: &[u8]) -> bool {
    return *TI.get(r).unwrap().get(a).unwrap()
}

fn get_info_float(rec: &rust_htslib::bcf::Record, field: &str) -> f32 {
    return rec.info(field.as_bytes()).float().unwrap().unwrap()[0]
}

fn is_snv(gt: i32, alleles: &Vec<&[u8]>) -> bool {
    return alleles[0].len() == alleles[gt as usize].len() && alleles[gt as usize].len() == 1
}

fn is_hom(gt: &[i32]) -> bool {
    return gt[0] == gt[1]
}

pub fn qc(
    input: Option<&str>,
    output: Option<&str>,
    field: Option<&str>,
    rare: &f32,
    gq: &i32,
    threads: &usize,
) {
    let mut bcf = bcfutils::get_rdr(input);
    bcf.set_threads(threads.clone())
        .expect("unable to set reader threads");

    let hdrv = bcf.header().clone();
    let n_samples = usize::try_from(hdrv.sample_count()).unwrap();
    let samples = hdrv.samples();

    let mut metrics = Metrics::new(n_samples);
    let mut rare_metrics = Metrics::new(n_samples);
    let mut qual_metrics = Metrics::new(n_samples);
    let mut rare_qual_metrics = Metrics::new(n_samples);

    let ff = field.expect("need to specify INFO field to get population frequency from");

    for record_result in bcf.records() {
        let record = record_result.expect("fail to read record");
        let pf = get_info_float(&record, ff);
        let alleles = record.alleles();

        let gqs = record.format(b"GQ").integer().expect("Couldn't retrieve GQ field");
        let gts = record.format(b"GT").integer().expect("Couldn't retrieve GT field");

        
        for sidx in 0..n_samples {

            //update if variant
            if gts[sidx][1] >= 4 {
                metrics.update(sidx, gts[sidx], &alleles);
                if pf <= *rare {
                    rare_metrics.update(sidx, gts[sidx], &alleles);
                }
                if gqs[sidx][0] >= *gq {
                    qual_metrics.update(sidx, gts[sidx], &alleles);
                }
                if pf <= *rare && gqs[sidx][0] >= *gq {
                    rare_qual_metrics.update(sidx, gts[sidx], &alleles)
                }
            }
        }
    }
    
    let iowtr: Box<dyn io::Write> = match output {
        None => Box::new(io::stdout()),
        Some("-") => Box::new(io::stdout()),
        Some(o) => Box::new(File::create(o).expect("unable to create output file"))
    };
    let mut wtr = csv::Writer::from_writer(iowtr);
    wtr.write_record(&[
        "sample", "n_var", "n_indels", "n_snvs", "n_hets", "n_hom_alts", "n_tis", "n_tvs", "ti_tv",
        format!("gq{gq}_n_var").as_str(), format!("gq{gq}_n_indels").as_str(), format!("gq{gq}_n_snvs").as_str(), format!("gq{gq}_n_hets").as_str(), format!("gq{gq}_n_hom_alts").as_str(), format!("gq{gq}_n_tis").as_str(), format!("gq{gq}_n_tvs").as_str(), format!("gq{gq}_ti_tv").as_str(),
        format!("{ff}{rare}_n_var").as_str(), format!("{ff}{rare}_n_indels").as_str(), format!("{ff}{rare}_n_snvs").as_str(), format!("{ff}{rare}_n_hets").as_str(), format!("{ff}{rare}_n_hom_alts").as_str(), format!("{ff}{rare}_n_tis").as_str(), format!("{ff}{rare}_n_tvs").as_str(), format!("{ff}{rare}_ti_tv").as_str(),
        format!("{ff}{rare}_gq{gq}_n_var").as_str(), format!("{ff}{rare}_gq{gq}_n_indels").as_str(), format!("{ff}{rare}_gq{gq}_n_snvs").as_str(), format!("{ff}{rare}_gq{gq}_n_hets").as_str(), format!("{ff}{rare}_gq{gq}_n_hom_alts").as_str(), format!("{ff}{rare}_gq{gq}_n_tis").as_str(), format!("{ff}{rare}_gq{gq}_n_tvs").as_str(), format!("{ff}{rare}_gq{gq}_ti_tv").as_str(),
    ]).expect("unable to write head row");
    for sidx in 0..n_samples {
        let ti_tv = metrics.n_tis[sidx] as f32 / metrics.n_tvs[sidx] as f32;
        let qual_ti_tv = qual_metrics.n_tis[sidx] as f32 / qual_metrics.n_tvs[sidx] as f32;
        let rare_ti_tv = rare_metrics.n_tis[sidx] as f32 / rare_metrics.n_tvs[sidx] as f32;
        let rare_qual_ti_tv = rare_qual_metrics.n_tis[sidx] as f32 / rare_qual_metrics.n_tvs[sidx] as f32;
        wtr.write_record(&[
            str::from_utf8(&samples[sidx]).unwrap(),
            &metrics.n_vars[sidx].to_string()[..],
            &metrics.n_indels[sidx].to_string()[..],
            &metrics.n_snvs[sidx].to_string()[..],
            &metrics.n_hets[sidx].to_string()[..],
            &metrics.n_hom_alts[sidx].to_string()[..],
            &metrics.n_tis[sidx].to_string()[..],
            &metrics.n_tvs[sidx].to_string()[..],
            &ti_tv.to_string()[..],
            &qual_metrics.n_vars[sidx].to_string()[..],
            &qual_metrics.n_indels[sidx].to_string()[..],
            &qual_metrics.n_snvs[sidx].to_string()[..],
            &qual_metrics.n_hets[sidx].to_string()[..],
            &qual_metrics.n_hom_alts[sidx].to_string()[..],
            &qual_metrics.n_tis[sidx].to_string()[..],
            &qual_metrics.n_tvs[sidx].to_string()[..],
            &qual_ti_tv.to_string()[..],
            &rare_metrics.n_vars[sidx].to_string()[..],
            &rare_metrics.n_indels[sidx].to_string()[..],
            &rare_metrics.n_snvs[sidx].to_string()[..],
            &rare_metrics.n_hets[sidx].to_string()[..],
            &rare_metrics.n_hom_alts[sidx].to_string()[..],
            &rare_metrics.n_tis[sidx].to_string()[..],
            &rare_metrics.n_tvs[sidx].to_string()[..],
            &rare_ti_tv.to_string()[..],
            &rare_qual_metrics.n_vars[sidx].to_string()[..],
            &rare_qual_metrics.n_indels[sidx].to_string()[..],
            &rare_qual_metrics.n_snvs[sidx].to_string()[..],
            &rare_qual_metrics.n_hets[sidx].to_string()[..],
            &rare_qual_metrics.n_hom_alts[sidx].to_string()[..],
            &rare_qual_metrics.n_tis[sidx].to_string()[..],
            &rare_qual_metrics.n_tvs[sidx].to_string()[..],
            &rare_qual_ti_tv.to_string()[..],
            ]).expect("unable to write row");
    }
}