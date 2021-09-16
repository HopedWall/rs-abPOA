use crate::abpoa::{
    abpoa_add_graph_edge, abpoa_add_graph_node, abpoa_align_sequence_to_graph, abpoa_dump_pog,
    abpoa_init, abpoa_init_para, abpoa_msa, abpoa_para_t, abpoa_post_set_para, abpoa_res_t,
    abpoa_t, free, strdup, ABPOA_CDEL, ABPOA_CDIFF, ABPOA_CHARD_CLIP, ABPOA_CINS, ABPOA_CMATCH,
    ABPOA_CSOFT_CLIP, ABPOA_SINK_NODE_ID, ABPOA_SRC_NODE_ID, FILE,
};
use rayon::prelude::*;
use std::collections::HashMap;
use std::ffi::CString;
use std::os::raw::{c_char, c_int};
use std::ptr;

pub struct AbpoaAligner {
    ab: *mut abpoa_t,
    abpt: *mut abpoa_para_t,

    // NOTE: the following only work when adding the nodes manually!!!
    // TODO: fix this
    n_nodes: usize,
    nodes: Vec<Vec<i32>>,
    // this does not consider the initial and final edge, however this should not
    // cause any issue
    edges: Vec<(usize, usize)>,
    edges_abpoa: Vec<(i32, i32)>,
}

pub struct AbpoaMSA {
    pub msa_length: usize,
    pub n_seqs: usize,
    // TODO: maybe this should be a map {seq_id : aln}?
    pub msa: Vec<String>,
}

impl AbpoaMSA {
    fn new() -> Self {
        AbpoaMSA {
            msa_length: 0,
            n_seqs: 0,
            msa: Vec::new(),
        }
    }

    fn new_from_alignment(msa: Vec<String>, n_seqs: usize, msa_length: usize) -> Self {
        AbpoaMSA {
            msa_length,
            n_seqs,
            msa,
        }
    }
}

pub struct AbpoaCons {
    pub cons_length: usize,
    pub cons: String,
}

impl AbpoaCons {
    fn new() -> Self {
        AbpoaCons {
            cons_length: 0,
            cons: String::new(),
        }
    }

    fn new_from_cons(cons: String) -> Self {
        AbpoaCons {
            cons_length: cons.len(),
            cons,
        }
    }
}

#[derive(Debug)]
pub struct AbpoaAlignmentResult {
    pub cigar: String,
    pub abpoa_nodes: Vec<u64>,
    pub graph_nodes: Vec<usize>,
    pub node_s: i32,
    pub node_e: i32,
    pub query_s: i32,
    pub query_e: i32,
    pub n_aligned_bases: i32,
    pub n_matched_bases: i32,
    pub best_score: i32,
    pub cigar_vec: Vec<char>,
}

impl AbpoaAlignmentResult {
    pub fn new() -> Self {
        AbpoaAlignmentResult {
            cigar: "".to_string(),
            abpoa_nodes: vec![],
            graph_nodes: vec![],
            node_s: 0,
            node_e: 0,
            query_s: 0,
            query_e: 0,
            n_aligned_bases: 0,
            n_matched_bases: 0,
            best_score: 0,
            cigar_vec: vec![],
        }
    }

    pub fn new_with_params(
        cigar: &str,
        abpoa_nodes: Vec<u64>,
        graph_nodes: Vec<usize>,
        res: &abpoa_res_t,
        cigar_vec: Vec<char>,
    ) -> Self {
        AbpoaAlignmentResult {
            cigar: cigar.to_string(),
            abpoa_nodes,
            graph_nodes,
            node_s: res.node_s,
            node_e: res.node_e,
            query_s: res.query_s,
            query_e: res.query_e,
            n_aligned_bases: res.n_aln_bases,
            n_matched_bases: res.n_matched_bases,
            best_score: res.best_score,
            cigar_vec,
        }
    }
}

impl AbpoaAligner {
    pub unsafe fn new() -> Self {
        AbpoaAligner {
            ab: abpoa_init(),
            abpt: abpoa_init_para(),
            n_nodes: 0,
            nodes: vec![],
            edges: vec![],
            edges_abpoa: vec![],
        }
    }

    // Initializes the aligner with the example.c params,
    // this is useful for debugging against the example.c
    // file in the original abpoa's repo
    pub unsafe fn new_with_example_params() -> Self {
        let mut aligner = AbpoaAligner::new();

        aligner.set_out_msa(true);
        aligner.set_out_cons(true);
        aligner.set_w(6);
        aligner.set_k(9);
        aligner.set_min_w(10);
        aligner.set_progressive_poa(true);
        aligner.set_post_para();

        aligner
    }

    pub unsafe fn set_out_msa(&mut self, val: bool) {
        (*self.abpt).set_out_msa(val as u8);
    }

    pub unsafe fn set_out_cons(&mut self, val: bool) {
        (*self.abpt).set_out_cons(val as u8);
    }

    pub unsafe fn set_progressive_poa(&mut self, val: bool) {
        (*self.abpt).set_progressive_poa(val as u8);
    }

    pub unsafe fn set_w(&mut self, w: u8) {
        (*self.abpt).w = w as c_int;
    }

    pub unsafe fn set_k(&mut self, k: u8) {
        (*self.abpt).k = k as c_int;
    }

    pub unsafe fn set_min_w(&mut self, min_w: u8) {
        (*self.abpt).min_w = min_w as c_int;
    }

    pub unsafe fn set_post_para(&mut self) {
        abpoa_post_set_para(self.abpt);
    }

    pub unsafe fn reset_aligner(&mut self) {
        (*(*self.ab).abs).n_seq = 0;
    }

    // NOTE: Rust does not support static fields, using const is the closest thing to that
    // see: https://stackoverflow.com/a/48972982
    const NT4_TABLE: [u8; 256] = [
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, /*'-'*/
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    ];
    const ALN_ALPHABET: [char; 6] = ['A', 'C', 'G', 'T', 'N', '-'];
    const CONS_ALPHABET: [char; 5] = ['A', 'C', 'G', 'T', 'N'];

    pub fn convert_seq_to_bseq(seq: &str) -> Vec<u8> {
        seq.chars()
            .map(|c| *(AbpoaAligner::NT4_TABLE).get(c as usize).unwrap())
            .collect()
    }

    pub unsafe fn align_seqs(&self, seqs: &Vec<&str>) -> AbpoaMSA {
        // Get the number of input sequences
        let n_seqs: c_int = seqs.len() as c_int;

        // Create a Vec with the sequences' length
        let mut seq_lens: Vec<c_int> = seqs.iter().map(|s| s.len() as c_int).collect();

        // Generate bseqs
        let mut bseqs_val: Vec<Vec<u8>> = seqs
            .into_iter()
            .map(|s| AbpoaAligner::convert_seq_to_bseq(s))
            .collect();

        let mut bseqs: Vec<*mut u8> = bseqs_val.iter_mut().map(|s| s.as_mut_ptr()).collect();

        // Now perform the alignment
        let mut cons_seq: *mut *mut u8 = ptr::null_mut();
        let mut cons_c: *mut *mut c_int = ptr::null_mut();
        let mut cons_l: *mut c_int = ptr::null_mut();
        let mut cons_n: c_int = 0;
        let mut msa_seq: *mut *mut u8 = ptr::null_mut();
        let mut msa_l: c_int = 0;
        let out: *mut FILE = ptr::null_mut(); //stdout;

        abpoa_msa(
            self.ab,
            self.abpt,
            n_seqs,
            ptr::null_mut(),
            seq_lens.as_mut_ptr(),
            bseqs.as_mut_ptr(),
            out,
            &mut cons_seq,
            &mut cons_c,
            &mut cons_l,
            &mut cons_n,
            &mut msa_seq,
            &mut msa_l,
        );

        // Read the alignment's results
        let mut msa: Vec<String> = Vec::new();
        for i in 0..n_seqs {
            let mut curr_aln = String::with_capacity(msa_l as usize);
            let outer_pointer = *msa_seq.add((i) as usize);
            for j in 0..msa_l {
                let inner_pointer = *(outer_pointer.add(j as usize));
                curr_aln.push(
                    *AbpoaAligner::ALN_ALPHABET
                        .get(inner_pointer as usize)
                        .unwrap(),
                );
            }
            msa.push(curr_aln);
        }

        AbpoaMSA::new_from_alignment(msa, n_seqs as usize, msa_l as usize)
    }

    pub unsafe fn consensus_from_seqs(&self, seqs: &Vec<&str>) -> AbpoaCons {
        // Get the number of input sequences
        let n_seqs: c_int = seqs.len() as c_int;

        // Create a Vec with the sequences' length
        let mut seq_lens: Vec<c_int> = seqs.iter().map(|s| s.len() as c_int).collect();

        // Generate bseqs
        let mut bseqs_val: Vec<Vec<u8>> = seqs
            .into_iter()
            .map(|s| AbpoaAligner::convert_seq_to_bseq(s))
            .collect();

        let mut bseqs: Vec<*mut u8> = bseqs_val.iter_mut().map(|s| s.as_mut_ptr()).collect();

        // Now perform the alignment
        let mut cons_seq: *mut *mut u8 = ptr::null_mut();
        let mut cons_c: *mut *mut c_int = ptr::null_mut();
        let mut cons_l: *mut c_int = ptr::null_mut();
        let mut cons_n: c_int = 0;
        let mut msa_seq: *mut *mut u8 = ptr::null_mut();
        let mut msa_l: c_int = 0;
        let out: *mut FILE = ptr::null_mut(); //stdout;

        abpoa_msa(
            self.ab,
            self.abpt,
            n_seqs,
            ptr::null_mut(),
            seq_lens.as_mut_ptr(),
            bseqs.as_mut_ptr(),
            out,
            &mut cons_seq,
            &mut cons_c,
            &mut cons_l,
            &mut cons_n,
            &mut msa_seq,
            &mut msa_l,
        );

        // Read the consensus
        let mut cons = String::with_capacity(*cons_l as usize);
        for i in 0..cons_n {
            let offset = *cons_l.add((i) as usize);
            for j in 0..offset {
                let outer_pointer = *cons_seq.add(i as usize);
                let inner_pointer = *(outer_pointer.add(j as usize));
                cons.push(
                    *AbpoaAligner::CONS_ALPHABET
                        .get(inner_pointer as usize)
                        .unwrap(),
                );
            }
        }

        AbpoaCons::new_from_cons(cons)
    }

    pub unsafe fn print_aln_to_dot(&mut self, path: &str) {
        // Build a C String to store path
        let c_str = CString::new(path).unwrap();

        (*self.abpt).out_pog = strdup(c_str.as_ptr() as *const c_char);
        abpoa_dump_pog(self.ab, self.abpt);
    }

    pub unsafe fn add_nodes_from_seq(&mut self, seq: &str) {
        let bseq: Vec<u8> = seq
            .chars()
            .map(|c| *(AbpoaAligner::NT4_TABLE).get(c as usize).unwrap())
            .collect();

        // First add the nodes to the graph
        // NOTE: in abpoa, each node has length 1 (i.e. a single nucleotide)
        let ids: Vec<i32> = bseq
            .into_iter()
            .map(|s| abpoa_add_graph_node((*self.ab).abg, s))
            .collect();

        //Then add the edges between said nodes
        ids.windows(2).for_each(|w| {
            abpoa_add_graph_edge(
                (*self.ab).abg,
                *w.get(0).unwrap(),
                *w.get(1).unwrap(),
                0,
                1,
                0,
                0,
                0,
            );
            self.edges_abpoa
                .push((*w.get(0).unwrap(), *w.get(1).unwrap()));
        });

        // Update wrapper data
        self.n_nodes += seq.len();
        self.nodes.push(ids.clone());
    }

    pub unsafe fn add_node(&mut self, base: u8) -> i32 {
        abpoa_add_graph_node((*self.ab).abg, base)
    }

    pub unsafe fn add_edge(&mut self, from_node_id: i32, to_node_id: i32) {
        abpoa_add_graph_edge((*self.ab).abg, from_node_id, to_node_id, 0, 1, 0, 0, 0);
        self.edges_abpoa.push((from_node_id, to_node_id));
    }

    pub unsafe fn add_nodes_edges(&mut self, nodes: &Vec<&str>, edges: &Vec<(usize, usize)>) {
        // Add nodes
        nodes.iter().for_each(|n| self.add_nodes_from_seq(n));

        // Add edges between nodes
        edges.iter().for_each(|e| {
            self.edges.push((e.0, e.1));
            let last_of_start_node = self.nodes.get(e.0).unwrap().last().unwrap().clone();
            let first_of_end_node = self.nodes.get(e.1).unwrap().first().unwrap().clone();
            self.add_edge(last_of_start_node, first_of_end_node);
        });

        if self.n_nodes > 0 {
            let heads: Vec<i32> = self.find_heads();
            //println!("Heads are: {:#?}", heads);
            for head in heads {
                // Add initial edge -- ABPOA_SRC_NODE_ID has node id 0
                abpoa_add_graph_edge(
                    (*self.ab).abg,
                    ABPOA_SRC_NODE_ID as i32,
                    head,
                    0,
                    1,
                    0,
                    0,
                    0,
                );
                self.edges_abpoa.push((ABPOA_SRC_NODE_ID as i32, head));
            }

            let tails: Vec<i32> = self.find_tails();
            //println!("Tails are: {:#?}", tails);
            for tail in tails {
                // Add initial edge -- ABPOA_SRC_NODE_ID has node id 0
                abpoa_add_graph_edge(
                    (*self.ab).abg,
                    tail,
                    ABPOA_SINK_NODE_ID as i32,
                    0,
                    1,
                    0,
                    0,
                    0,
                );
                self.edges_abpoa.push((tail, ABPOA_SINK_NODE_ID as i32));
            }
        }
    }

    pub unsafe fn find_heads(&self) -> Vec<i32> {
        let heads = self
            .nodes
            .clone()
            .into_iter()
            .flatten()
            .filter(|abpoa_node| self.is_head(abpoa_node))
            .collect();

        heads
    }

    fn is_head(&self, node: &i32) -> bool {
        let node_appears_as_end_edge = self.edges_abpoa.iter().any(|x| x.1 == *node);
        return !node_appears_as_end_edge;
    }

    pub unsafe fn find_tails(&self) -> Vec<i32> {
        let tails = self
            .nodes
            .clone()
            .into_iter()
            .flatten()
            .filter(|abpoa_node| self.is_tail(abpoa_node))
            .collect();

        tails
    }

    fn is_tail(&self, node: &i32) -> bool {
        let node_appears_as_start_edge = self.edges_abpoa.iter().any(|x| x.0 == *node);
        return !node_appears_as_start_edge;
    }

    pub unsafe fn align_sequence(&mut self, seq: &str) -> AbpoaAlignmentResult {
        let mut bseq: Vec<u8> = AbpoaAligner::convert_seq_to_bseq(seq);

        let mut res = abpoa_res_t {
            n_cigar: 0,
            m_cigar: 0,
            graph_cigar: vec![].as_mut_ptr(),
            node_s: 0,
            node_e: 0,
            query_s: 0,
            query_e: 0,
            n_aln_bases: 0,
            n_matched_bases: 0,
            best_score: 0,
        };

        abpoa_align_sequence_to_graph(
            self.ab,
            self.abpt,
            bseq.as_mut_ptr(),
            seq.len() as i32,
            &mut res,
        );

        // Build a HashMap having as keys the abpoa_ids and as values
        // the abstraction_id, this will be useful when doing the conversion
        // TODO: would like to use rayon but it does not like the .get()
        let mut abpoa_id_to_abstraction_id: HashMap<u64, usize> = HashMap::new();
        for abstraction_id in 0..self.nodes.len() {
            for abpoa_id in self.nodes.get(abstraction_id).unwrap() {
                abpoa_id_to_abstraction_id.insert(*abpoa_id as u64, abstraction_id);
            }
        }

        // Create variables to store aln result
        let mut abpoa_ids: Vec<u64> = Vec::new();
        let mut cigar_vec: Vec<char> = Vec::new();

        // Navigate the cigar
        for i in 0..res.n_cigar {
            let curr_cigar = res.graph_cigar.add(i as usize);
            let op = (*curr_cigar & 0xf) as u32;

            let op_char = match op {
                ABPOA_CMATCH => 'M',
                ABPOA_CINS => 'I',
                ABPOA_CDEL => 'D',
                ABPOA_CDIFF => 'X',
                ABPOA_CSOFT_CLIP => 'S',
                ABPOA_CHARD_CLIP => 'H',
                _ => ' ',
            };

            let node_id = (*curr_cigar >> 34) & 0x3fffffff;

            // Necessary because sometimes abpoa returns weird nodes
            // TODO: figure out why this happens
            if abpoa_id_to_abstraction_id.contains_key(&node_id) && op_char != ' ' {
                abpoa_ids.push(node_id);
                cigar_vec.push(op_char);
            }
        }

        // Convert abpoa_ids to abstraction_ids
        let graph_ids: Vec<usize> = abpoa_ids
            .par_iter()
            .filter_map(|id| abpoa_id_to_abstraction_id.get(id))
            .map(|id| *id)
            .collect();

        // Compact cigar
        let mut cigar_string = String::new();
        if !cigar_vec.is_empty() {
            let mut last_char = ' ';
            let mut count = 0;

            for char in &cigar_vec {
                if *char != last_char {
                    if last_char != ' ' {
                        // This is the initial delimiter for this loop
                        cigar_string.push_str(&mut format!("{}{}", count, last_char));
                    }
                    last_char = *char;
                    count = 1;
                } else {
                    count += 1;
                }
            }

            cigar_string.push_str(&mut format!("{}{}", count, last_char));
        }

        assert_eq!(abpoa_ids.len(), graph_ids.len());

        AbpoaAlignmentResult::new_with_params(
            cigar_string.as_str(),
            abpoa_ids,
            graph_ids,
            &res,
            cigar_vec,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::abpoa::{abpoa_generate_gfa, stdout};

    // ----- Check msa and consensus ("black-box" version) -----
    #[test]
    fn test_aln() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();

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
                "CGTCAATCTAATCGAAGCATACGCGGCAGAGCCGTCTACCTCGGCAATCACGT",
            ]
            .to_vec();

            let aln = aligner.align_seqs(&seqs);

            //aligner.print_aln_to_dot("example.png");

            //aligner.reset_aligner();

            //println!("MSA: {:#?}", aln.msa);
            assert_eq!(aln.n_seqs, seqs.len());
            assert_eq!(aln.msa_length, 75);
        }
    }

    #[test]
    fn test_cons() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();

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
                "CGTCAATCTAATCGAAGCATACGCGGCAGAGCCGTCTACCTCGGCAATCACGT",
            ]
            .to_vec();

            let cons = aligner.consensus_from_seqs(&seqs);
            //aligner.print_aln_to_dot("example.png");

            //aligner.reset_aligner();
            //assert_eq!(aln.n_seqs, seqs.len());
            //assert_eq!(aln.msa_length, 75);
        }
    }

    // ----- Test basic abstraction functionalities -----
    #[test]
    fn test_add_nodes() {
        unsafe {
            let mut aligner = AbpoaAligner::new();
            aligner.add_nodes_from_seq("ACGT");
            assert_eq!(aligner.nodes.first().unwrap().len(), 4);
            assert_eq!(aligner.n_nodes, 4);
        }
    }

    #[test]
    fn test_add_nodes_and_edges() {
        unsafe {
            let mut aligner = AbpoaAligner::new();
            aligner.add_nodes_edges(&vec!["ACG", "GCT", "TAT"], &vec![(0, 1), (0, 2)]);
            assert_eq!(aligner.nodes.len(), 3);
            assert_eq!(aligner.n_nodes, 9);
            assert_eq!(aligner.edges.len(), 2);
        }
    }

    // ----- Test graph-to-seq alignment -----
    #[test]
    fn test_alignment_1() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&vec!["ACT"], &vec![]);
            let res = aligner.align_sequence("ACT");

            assert_eq!(res.cigar, String::from("3M"));
            // These are the nodes in abpoa (remember, nodes are 1-base only!)
            assert_eq!(res.abpoa_nodes, vec![2, 3, 4]);
            // These are the nodes in our graph abstraction (nodes can have length > 1)
            assert_eq!(res.graph_nodes, vec![0, 0, 0]);
        }
    }

    #[test]
    fn test_alignment_2() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&vec!["ACT"], &vec![]);
            let res = aligner.align_sequence("T");

            assert_eq!(res.cigar, String::from("1M"));
            // These are the nodes in abpoa (remember, nodes are 1-base only!)
            assert_eq!(res.abpoa_nodes, vec![4]);
            // These are the nodes in our graph abstraction (nodes can have length > 1)
            assert_eq!(res.graph_nodes, vec![0]);
        }
    }

    #[test]
    fn test_alignment_3() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&vec!["ACG"], &vec![]);
            let res = aligner.align_sequence("ATG");

            //Unexpected but same behavior in C ver
            assert_eq!(res.cigar, String::from("3M"));
            // These are the nodes in abpoa (remember, nodes are 1-base only!)
            assert_eq!(res.abpoa_nodes, vec![2, 3, 4]);
            // These are the nodes in our graph abstraction (nodes can have length > 1)
            assert_eq!(res.graph_nodes, vec![0, 0, 0]);
        }
    }

    #[test]
    fn test_alignment_3_manual() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&vec!["ACG"], &vec![]);
            let res = aligner.align_sequence("ATG");

            // Check against the "manual" version (= less wrapper abstractions used, closer
            // to the original C impl.)
            let res_manual = manual_test_single_node("ACG", "ATG");
            assert_eq!(res.cigar, res_manual.cigar);
            assert_eq!(res.abpoa_nodes, res_manual.abpoa_nodes);

            // Makes no sense to compare res.graph_nodes because it is an abstraction
            // only available in Rust
        }
    }

    #[test]
    fn test_alignment_4_multiple() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&vec!["ACG", "AAA"], &vec![(0, 1)]);
            let res = aligner.align_sequence("ATG");

            //Unexpected but same behavior in C ver
            assert_eq!(res.cigar, String::from("3M3D"));
            // These are the nodes in abpoa (remember, nodes are 1-base only!)
            assert_eq!(res.abpoa_nodes, vec![2, 3, 4, 5, 6, 7]);
            // These are the nodes in our graph abstraction (nodes can have length > 1)
            assert_eq!(res.graph_nodes, vec![0, 0, 0, 1, 1, 1]);
        }
    }

    #[test]
    fn test_alignment_4_multiple_manual() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&vec!["ACG", "AAA"], &vec![(0, 1)]);
            let res = aligner.align_sequence("ATG");

            // Check against the "manual" version (= less wrapper abstractions used, closer
            // to the original C impl.)
            let res_manual = manual_test_multiple_nodes("ACG", "AAA", "ATG");
            assert_eq!(res.cigar, res_manual.cigar);
            assert_eq!(res.abpoa_nodes, res_manual.abpoa_nodes);

            // Makes no sense to compare res.graph_nodes because it is an abstraction
            // only available in Rust
        }
    }

    fn manual_test_single_node(seq: &str, query: &str) -> AbpoaAlignmentResult {
        unsafe {
            let mut aligner = AbpoaAligner::new();

            aligner.set_out_cons(true);
            aligner.set_out_cons(true);
            aligner.set_w(6);
            aligner.set_k(9);
            aligner.set_min_w(10);
            aligner.set_progressive_poa(true);

            abpoa_post_set_para(aligner.abpt);

            let bseq: Vec<u8> = seq
                .chars()
                .map(|c| *(AbpoaAligner::NT4_TABLE).get(c as usize).unwrap())
                .collect();

            let mut ids: Vec<i32> = Vec::new();

            for b in bseq {
                let c_id = aligner.add_node(b);
                ids.push(c_id);
            }
            //println!("ids: {:#?}", ids);

            let mut prev_node_id = ABPOA_SRC_NODE_ID as i32;
            let mut curr_node_id: i32 = 0;

            for i in 0..ids.len() {
                curr_node_id = *ids.get(i).unwrap();
                aligner.add_edge(prev_node_id, curr_node_id);
                prev_node_id = curr_node_id;
            }
            aligner.add_edge(curr_node_id, ABPOA_SINK_NODE_ID as i32);

            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);

            let mut query_bseq: Vec<u8> = query
                .chars()
                .map(|c| *(AbpoaAligner::NT4_TABLE).get(c as usize).unwrap())
                .collect();

            let mut res = abpoa_res_t {
                n_cigar: 0,
                m_cigar: 0,
                graph_cigar: vec![].as_mut_ptr(),
                node_s: 0,
                node_e: 0,
                query_s: 0,
                query_e: 0,
                n_aln_bases: 0,
                n_matched_bases: 0,
                best_score: 0,
            };
            abpoa_align_sequence_to_graph(
                aligner.ab,
                aligner.abpt,
                query_bseq.as_mut_ptr(),
                query.len() as i32,
                &mut res,
            );

            // Create variables to store aln result
            let mut abpoa_ids: Vec<u64> = Vec::new();
            let mut graph_ids: Vec<usize> = Vec::new();
            let mut cigar_vec: Vec<char> = Vec::new();

            // Navigate the cigar
            let mut op: u32 = 0;
            let mut op_char = ' ';
            let mut node_id: u64 = 0;
            for i in 0..res.n_cigar {
                let curr_cigar = res.graph_cigar.add(i as usize);
                op = (*curr_cigar & 0xf) as u32;

                op_char = match op {
                    ABPOA_CMATCH => 'M',
                    ABPOA_CINS => 'I',
                    ABPOA_CDEL => 'D',
                    ABPOA_CDIFF => 'X',
                    ABPOA_CSOFT_CLIP => 'S',
                    ABPOA_CHARD_CLIP => 'H',
                    _ => ' ',
                };

                node_id = (*curr_cigar >> 34) & 0x3fffffff;
                abpoa_ids.push(node_id);
                cigar_vec.push(op_char);
            }

            //println!("Ids {:#?}", abpoa_ids);
            //println!("Cigar vec {:#?}", cigar_vec);
            // Compact cigar
            let mut cigar_string = String::new();
            if !cigar_vec.is_empty() {
                let mut last_char = ' ';
                let mut count = 0;

                for char in &cigar_vec {
                    if *char != last_char {
                        if last_char != ' ' {
                            // This is the initial delimiter
                            cigar_string.push_str(&mut format!("{}{}", count, last_char));
                        }
                        last_char = *char;
                        count = 1;
                    } else {
                        count += 1;
                    }
                }

                cigar_string.push_str(&mut format!("{}{}", count, last_char));
            }

            AbpoaAlignmentResult::new_with_params(
                cigar_string.as_str(),
                abpoa_ids,
                graph_ids,
                &res,
                cigar_vec,
            )
        }
    }

    fn manual_test_multiple_nodes(seq: &str, seq2: &str, query: &str) -> AbpoaAlignmentResult {
        unsafe {
            let mut aligner = AbpoaAligner::new();

            aligner.set_out_cons(true);
            aligner.set_out_cons(true);
            aligner.set_w(6);
            aligner.set_k(9);
            aligner.set_min_w(10);
            aligner.set_progressive_poa(true);

            abpoa_post_set_para(aligner.abpt);

            let bseq: Vec<u8> = seq
                .chars()
                .map(|c| *(AbpoaAligner::NT4_TABLE).get(c as usize).unwrap())
                .collect();

            let mut ids: Vec<i32> = Vec::new();

            for b in bseq {
                let c_id = aligner.add_node(b);
                ids.push(c_id);
            }
            //println!("ids: {:#?}", ids);

            let mut prev_node_id = ABPOA_SRC_NODE_ID as i32;
            let mut curr_node_id: i32 = 0;

            for i in 0..ids.len() {
                curr_node_id = *ids.get(i).unwrap();
                aligner.add_edge(prev_node_id, curr_node_id);
                prev_node_id = curr_node_id;
            }

            // Add node 2
            let bseq2: Vec<u8> = seq2
                .chars()
                .map(|c| *(AbpoaAligner::NT4_TABLE).get(c as usize).unwrap())
                .collect();

            let mut ids2: Vec<i32> = Vec::new();
            for b in bseq2 {
                let c_id = aligner.add_node(b);
                ids2.push(c_id);
            }
            for i in 0..ids2.len() {
                curr_node_id = *ids2.get(i).unwrap();
                aligner.add_edge(prev_node_id, curr_node_id);
                prev_node_id = curr_node_id;
            }

            // Add final edge
            aligner.add_edge(curr_node_id, ABPOA_SINK_NODE_ID as i32);

            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);

            let mut query_bseq: Vec<u8> = query
                .chars()
                .map(|c| *(AbpoaAligner::NT4_TABLE).get(c as usize).unwrap())
                .collect();

            let mut res = abpoa_res_t {
                n_cigar: 0,
                m_cigar: 0,
                graph_cigar: vec![].as_mut_ptr(),
                node_s: 0,
                node_e: 0,
                query_s: 0,
                query_e: 0,
                n_aln_bases: 0,
                n_matched_bases: 0,
                best_score: 0,
            };
            abpoa_align_sequence_to_graph(
                aligner.ab,
                aligner.abpt,
                query_bseq.as_mut_ptr(),
                query.len() as i32,
                &mut res,
            );

            // Create variables to store aln result
            let mut abpoa_ids: Vec<u64> = Vec::new();
            let mut graph_ids: Vec<usize> = Vec::new();
            let mut cigar_vec: Vec<char> = Vec::new();

            // Navigate the cigar
            let mut op: u32 = 0;
            let mut op_char = ' ';
            let mut node_id: u64 = 0;
            for i in 0..res.n_cigar {
                let curr_cigar = res.graph_cigar.add(i as usize);
                op = (*curr_cigar & 0xf) as u32;

                op_char = match op {
                    ABPOA_CMATCH => 'M',
                    ABPOA_CINS => 'I',
                    ABPOA_CDEL => 'D',
                    ABPOA_CDIFF => 'X',
                    ABPOA_CSOFT_CLIP => 'S',
                    ABPOA_CHARD_CLIP => 'H',
                    _ => ' ',
                };

                node_id = (*curr_cigar >> 34) & 0x3fffffff;

                //println!("Node id {} type {}", node_id, op_char);

                abpoa_ids.push(node_id);
                cigar_vec.push(op_char);
            }

            // Compact cigar
            let mut cigar_string = String::new();
            if !cigar_vec.is_empty() {
                let mut last_char = ' ';
                let mut count = 0;

                for char in &cigar_vec {
                    if *char != last_char {
                        if last_char != ' ' {
                            // This is the initial delimiter
                            cigar_string.push_str(&mut format!("{}{}", count, last_char));
                        }
                        last_char = *char;
                        count = 1;
                    } else {
                        count += 1;
                    }
                }

                cigar_string.push_str(&mut format!("{}{}", count, last_char));
            }

            AbpoaAlignmentResult::new_with_params(
                cigar_string.as_str(),
                abpoa_ids,
                graph_ids,
                &res,
                cigar_vec,
            )
        }
    }

    #[test]
    fn add_complex_graph_1() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&vec!["ACG", "AAA"], &vec![(0, 1)]);
            abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
        }
    }

    #[test]
    fn add_complex_graph_2() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&vec!["ACG", "AAA", "CC"], &vec![(0, 1)]);
            abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
        }
    }

    #[test]
    fn add_complex_graph_3() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();
            let nodes: Vec<&str> = vec![
                "A", "G", "AAAT", "AA", "TTTCT", "GG", "AGTTCTAT", "A", "T", "ATAT", "A", "T",
            ];

            let edges: Vec<(usize, usize)> = vec![
                (0, 2),
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (4, 5),
                (4, 6),
                (5, 6),
                (6, 7),
                (6, 8),
                (7, 9),
                (8, 9),
                (9, 10),
                (9, 11),
            ];
            //(0,1), (7,8), (10,11)];
            aligner.add_nodes_edges(&nodes, &edges);
            //println!("Nodes: {:?}", aligner.nodes);
            //println!("Edges: {:?}", aligner.edges);
            //println!("Edges_abpoa: {:?}", aligner.edges_abpoa);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);

            let result = aligner.align_sequence("AAATTTGGCAT");
            //println!("Result is: {:#?}", result);

            assert_eq!(
                result.abpoa_nodes,
                vec![
                    2, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                    27, 28, 29, 30, 31
                ]
            );

            assert_eq!(
                result.graph_nodes,
                vec![0, 2, 2, 2, 2, 4, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 9, 9, 9, 9, 10]
            );
        }
    }
}
