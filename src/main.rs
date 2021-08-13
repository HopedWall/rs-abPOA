use ab_poa::abpoa::*;
use std::ptr::null;

fn main() {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{mem, io, ptr};
    use std::os::raw::{c_int, c_char};
    use std::alloc::alloc;

    #[test]
    fn test() {
        unsafe {
            let ab_poa = abpoa_init();
            let mut abpt = abpoa_init_para();

            (*abpt).set_out_msa(1);
            (*abpt).set_out_cons(1);
            (*abpt).w = 6;
            (*abpt).k = 9;
            (*abpt).min_w = 10;
            (*abpt).set_progressive_poa(1);

            abpoa_post_set_para(abpt);

            let seqs: [&str; 10] = [
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
            ];

            let _char256_table : [char; 256] = [
                'A', 'C', 'G', 'T',  'N', 'B', 'D', 'E',  'F', 'H', 'I',  'J', 'K', 'L', 'M', 'O',
                'P', 'Q', 'R', 'S',  'U', 'V', 'W', 'X',  'Y', 'Z', '*', '-',  '*', '*', '*', '*',
                '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',
                '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',
                '*', 'A', 'B', 'C',  'D', 'E', 'F', 'G',  'H', 'I', 'J', 'K',  'L', 'M', 'N', 'O',
                'P', 'Q', 'R', 'S',  'T', 'U', 'V', 'W',  'X', 'Y', 'Z', '*',  '*', '*', '*', '*',
                '*', 'A', 'B', 'C',  'D', 'E', 'F', 'G',  'H', 'I', 'J', 'K',  'L', 'M', 'N', 'O',
                'P', 'Q', 'R', 'S',  'T', 'U', 'V', 'W',  'X', 'Y', 'Z', '*',  '*', '*', '*', '*',
                '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',
                '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',
                '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',
                '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',
                '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',
                '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',
                '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',
                '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*'
            ];

             let _char26_table: [u8; 256] = [
                0,  1,  2,  3,   4,  5,  6,  7,   8,  9, 10, 11,  12, 13, 14, 15,
                16, 17, 18, 19,  20, 21, 22, 23,  24, 25, 26, 26,  26, 26, 26, 26,
                26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
                26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
                26,  0,  5,  1,   6,  7,  8,  2,   9, 10, 11, 12,  13, 14,  4, 15,
                16, 17, 18, 19,   3, 20, 21, 22,  23, 24, 25, 26,  26, 26, 26, 26,
                26,  0,  5,  1,   6,  7,  8,  2,   9, 10, 11, 12,  13, 14,  4, 15,
                16, 17, 18, 19,   3, 20, 21, 22,  23, 24, 25, 26,  26, 26, 26, 26,
                26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
                26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
                26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
                26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
                26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
                26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
                26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
                26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26
            ];


            // c_int is an alias for i32
            let n_seqs: c_int = 10;
            //let seq_lens: *mut c_int  = malloc((mem::size_of::<i32>() * n_seqs as usize) as u64) as *mut c_int;
            let mut seq_lens_val : [c_int; 10] = [0; 10];
            let seq_lens: *mut c_int  = seq_lens_val.as_mut_ptr();

            //let bseqs: *mut *mut u8 = malloc((mem::size_of::<u8>() * n_seqs as usize) as u64) as *mut *mut u8;
            let mut bseqs_val: [* mut u8; 10] = [ptr::null_mut();10];
            let bseqs = bseqs_val.as_mut_ptr();
            for i in 0..n_seqs {
                //seq_lens[i] = &strlen(seqs[i]);
                let mut curr_seq_len : &mut c_int = seq_lens_val.get_mut(i as usize).unwrap();
                *curr_seq_len = seqs.get(i as usize).unwrap().len() as c_int;
                //bseqs[i] = malloc(u8::size_of() * seq_lens[i]);
                let mut curr_numbers : Vec<u8> = Vec::new();
                for j in 0..*curr_seq_len {
                    curr_numbers.push(*_char26_table.get((i*16+j) as usize).unwrap());
                }
            }
            
            println!("Seq lens: {:#?}", seq_lens_val);

            let cons_seq: *mut *mut *mut u8 = ptr::null_mut();
            let cons_cov: *mut *mut *mut c_int = ptr::null_mut();
            let cons_l : *mut *mut c_int = ptr::null_mut();
            let cons_n: *mut c_int = ptr::null_mut();
            let msa_seq: *mut *mut *mut u8 = ptr::null_mut();
            let msa_l : *mut c_int = ptr::null_mut();

            //abpoa_msa(ab_poa, abpt, 10, ptr::null_mut(), seq_lens, bseqs, ptr::null_mut(), cons_seq, cons_cov, cons_l, cons_n, msa_seq, msa_l);

            //let msa_seq_val = msa_seq;
            /*
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
