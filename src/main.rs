use ab_poa::abpoa::*;
use std::ptr::null;

fn main() {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{mem, io, ptr};
    use std::os::raw::{c_int, c_char};
    use std::alloc::alloc;
    use std::ops::Deref;

    #[test]
    fn test() {
        unsafe {
            let ab = abpoa_init();
            let mut abpt = abpoa_init_para();

            (*abpt).set_out_msa(1);
            (*abpt).set_out_cons(1);
            (*abpt).w = 6;
            (*abpt).k = 9;
            (*abpt).min_w = 10;
            (*abpt).set_progressive_poa(1);

            //println!("abpt: {:#?}", *abpt);

            abpoa_post_set_para(abpt);

            let seqs: Vec<&str> = [
                "CGTCAATCTATCGAAGCATACGCGGGCAGAGCCGAAGACCTCGGCAATCCA",
                "CCACGTCAATCTATCGAAGCATACGCGGCAGCCGAACTCGACCTCGGCAATCAC",
                "CGTCAATCTATCGAAGCATACGCGGCAGAGCCCGGAAGACCTCGGCAATCAC",
                "CGTCAATGCTAGTCGAAGCAGCTGCGGCAGAGCCGAAGACCTCGGCAATCAC",
                "CGTCAATCTATCGAAGCATTCTACGCGGCAGAGCCGACCTCGGCAATCAC",
                "CGTCAATCTAGAAGCATACGCGGCAAGAGCCGAAGACCTCGGCCAATCAC",
                "CGTCAATCTATCGGTAAAGCATACGCTCTGTAGCCGAAGACCTCGGCAATCAC",
                "CGTCAATCTATCTTCAAGCATACGCGGCAGAGCCGAAGACCTCGGCAATC",
                "CGTCAATGGATCGAGTACGCGGCAGAGCCGAAGACCTCGGCAATCAC",
                "CGTCAATCTAATCGAAGCATACGCGGCAGAGCCGTCTACCTCGGCAATCACGT"
            ].to_vec();

            let _nt4_table : [u8; 256] = [
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
            ];


            // Get the number of input sequences
            let n_seqs: c_int = seqs.len() as c_int;

            // Create a Vec with the sequences' length
            let mut seq_lens_val : Vec<c_int> = vec![0;seqs.len()];
            let seq_lens: *mut c_int  = seq_lens_val.as_mut_ptr();

            // Create a matrix (bseqs_val) where:
            // - each row represents an input sequence (rows can have different lengths)
            // - each column represents a base as per the _nth4_table
            // Bseqs is a Vec of pointers each pointing to one row

            let mut bseqs_val_val : Vec<Vec<u8>> = Vec::new();
            for i in 0..n_seqs {
                bseqs_val_val.push(Vec::new());
            }
            let mut bseqs_val: Vec<*mut u8> = Vec::with_capacity(seqs.len());
            for i in 0..bseqs_val.len() {
                let mut bval = bseqs_val.get_mut(i).unwrap();
                *bval = bseqs_val_val.get_mut(i as usize).unwrap().as_mut_ptr();
            }

            let bseqs = bseqs_val.as_mut_ptr();
            for i in 0..n_seqs {
                let curr_seq : &str = seqs.get(i as usize).unwrap();
                let mut curr_seq_len : &mut c_int = seq_lens_val.get_mut(i as usize).unwrap();
                *curr_seq_len = curr_seq.len() as c_int;
                let mut curr_numbers : Vec<u8> = curr_seq.chars()
                                                .map(|c| *(_nt4_table).get(c as usize).unwrap())
                                                .collect();

                let mut bval : &mut Vec<u8> = bseqs_val_val.get_mut(i as usize).unwrap();
                *bval = curr_numbers;
            }

            // Test: print seq_lens using pointers
            for i in 0..n_seqs {
                println!("Seq_lens[{}] is: {}", i, *seq_lens.add(i as usize));
            }

            /*
            // Test: print bseqs using pointers
            for i in 0..n_seqs {
                let curr_bseq = bseqs;
                println!("value is: {}", **curr_bseq.add(i as usize));
            }
            for i in 0..n_seqs {
                let mut curr_bseq_val = bseqs_val.get_mut(i as usize).unwrap();

                let curr_seq_lens = *seq_lens_val.get(i as usize).unwrap();
                println!("Curr seq lens: {}", curr_seq_lens);
                for j in 0..curr_seq_lens {
                    println!("Bseq [{}][{}]: {:#?}", i, j, *(curr_bseq_val.add(j as usize)));
                }

                println!("");
            }
             */

            abpoa_msa(ab, abpt, n_seqs, ptr::null_mut(), seq_lens, bseqs, ptr::null_mut(),
                      ptr::null_mut(), ptr::null_mut(), ptr::null_mut(),
                      ptr::null_mut(), ptr::null_mut(), ptr::null_mut());

            /*
            let cons_seq: *mut *mut *mut u8 = ptr::null_mut();
            let cons_cov: *mut *mut *mut c_int = ptr::null_mut();
            let cons_l : *mut *mut c_int = ptr::null_mut();
            let cons_n: *mut c_int = &mut 0;
            let msa_seq: *mut *mut *mut u8 = ptr::null_mut();
            let msa_l : *mut c_int = &mut 0;

            abpoa_msa(ab_poa, abpt, 10, ptr::null_mut(), seq_lens, bseqs, ptr::null_mut(), cons_seq, cons_cov, cons_l, cons_n, msa_seq, msa_l);

            println!(">Multiple_sequence_alignment");
            for i in 0..n_seqs {
                print!("ACGTN-");
                for j in 0..*msa_l {
                    print!("{}", msa_seq[i][j])
                }
                println!();
            }
             */

        }
    }
}
