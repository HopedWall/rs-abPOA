use crate::abpoa::{abpoa_add_graph_edge, abpoa_add_graph_node, abpoa_align_sequence_to_graph, abpoa_dump_pog, abpoa_free, abpoa_free_para, abpoa_init, abpoa_init_para, abpoa_msa, abpoa_para_t, abpoa_post_set_para, abpoa_res_t, abpoa_t, free, strdup, ABPOA_CDEL, ABPOA_CDIFF, ABPOA_CHARD_CLIP, ABPOA_CINS, ABPOA_CMATCH, ABPOA_CSOFT_CLIP, ABPOA_SINK_NODE_ID, ABPOA_SRC_NODE_ID, FILE, ABPOA_LOCAL_MODE, ABPOA_GLOBAL_MODE};
//use rayon::prelude::*;
use std::collections::HashMap;
use std::ffi::{c_void, CString};
use std::os::raw::{c_char, c_int};
use std::ptr;

use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use log::{info, warn};
use std::env;
use std::iter::once;
use std::path::PathBuf;
use std::time::Instant;

pub struct AbpoaAligner {
    ab: *mut abpoa_t,
    abpt: *mut abpoa_para_t,

    // NOTE: the following only work when adding the nodes manually!!!
    // TODO: fix this
    n_nodes: usize,
    nodes: Vec<Vec<i32>>,
    nodes_str: HashMap<i32, char>,
    abpoa_id_to_abstraction_id: HashMap<i32, usize>,

    // this does not consider the initial and final edge, however this should not
    // cause any issue
    edges: Vec<(usize, usize)>,
    edges_abpoa: Vec<(i32, i32)>,
}

pub enum AbpoaAlignmentMode {
    Global,
    Local
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
    pub abpoa_nodes: Vec<i32>,
    pub node_s: i32,
    pub node_e: i32,
    pub graph_nodes: Vec<usize>,
    pub aln_start_offset: usize,
    pub aln_end_offset: usize,
    pub query_s: i32,
    pub query_e: i32,
    pub n_aligned_bases: i32,
    pub n_matched_bases: i32,
    pub best_score: i32,
    pub cigar_vec: Vec<char>,
    pub cs_string: String,
}

impl AbpoaAlignmentResult {
    pub fn new() -> Self {
        AbpoaAlignmentResult {
            cigar: "".to_string(),
            abpoa_nodes: vec![],
            node_s: 0,
            node_e: 0,
            graph_nodes: vec![],
            aln_start_offset: 0,
            aln_end_offset: 0,
            query_s: 0,
            query_e: 0,
            n_aligned_bases: 0,
            n_matched_bases: 0,
            best_score: 0,
            cigar_vec: vec![],
            cs_string: "".to_string(),
        }
    }

    pub fn new_with_params(
        cigar: &str,
        abpoa_nodes: Vec<i32>,
        graph_nodes: Vec<usize>,
        aln_start_offset: usize,
        aln_end_offset: usize,
        res: &abpoa_res_t,
        cigar_vec: Vec<char>,
        cs_string: &str,
    ) -> Self {
        AbpoaAlignmentResult {
            cigar: cigar.to_string(),
            abpoa_nodes,
            node_s: res.node_s,
            node_e: res.node_e,
            graph_nodes,
            aln_start_offset,
            aln_end_offset,
            query_s: res.query_s,
            query_e: res.query_e,
            n_aligned_bases: res.n_aln_bases,
            n_matched_bases: res.n_matched_bases,
            best_score: res.best_score,
            cigar_vec,
            cs_string: cs_string.to_string(),
        }
    }
}

impl AbpoaAligner {
    pub unsafe fn new() -> Self {

        let aln = AbpoaAligner {
            ab: abpoa_init(),
            abpt: abpoa_init_para(),
            n_nodes: 0,
            nodes: vec![],
            nodes_str: HashMap::new(),
            edges: vec![],
            edges_abpoa: vec![],
            abpoa_id_to_abstraction_id: HashMap::new(),
        };

        //println!("Init para: {:#?}", *aln.abpt);

        aln
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

    pub unsafe fn new_with_params(w: u8, k: u8, min_w: u8, abpoa_mode: AbpoaAlignmentMode) -> Self {
        let mut aligner = AbpoaAligner::new();

        aligner.set_w(w);
        aligner.set_k(k);
        aligner.set_min_w(min_w);

        match abpoa_mode {
            AbpoaAlignmentMode::Global => (*aligner.abpt).align_mode = ABPOA_GLOBAL_MODE as c_int,
            AbpoaAlignmentMode::Local => (*aligner.abpt).align_mode = ABPOA_LOCAL_MODE as c_int,
        };

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
    const NT256_TABLE: [char; 256] = [
        'A', 'C', 'G', 'T', 'N', '-', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', '-', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'A', 'N', 'C', 'N', 'N', 'N', 'G',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'T', 'T', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'A', 'N', 'C', 'N', 'N', 'N', 'G', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'T', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',
    ];
    const ALN_ALPHABET: [char; 6] = ['A', 'C', 'G', 'T', 'N', '-'];
    const CONS_ALPHABET: [char; 5] = ['A', 'C', 'G', 'T', 'N'];

    /*
    pub fn convert_seq_to_bseq(seq: &str) -> Vec<u8> {
        seq.chars()
            .map(|c| *(AbpoaAligner::NT4_TABLE).get(c as usize).unwrap())
            .collect()
    }
     */

    pub fn convert_seq_to_bseq(seq: &str) -> Vec<u8> {
        seq.chars()
            .map(|c| match c {
                'A' | 'a' => 0,
                'C' | 'c' => 1,
                'G' | 'g' => 2,
                'T' | 't' => 3,
                _ => panic!("Cannot add base to graph!"),
            })
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
        );

        // Read the alignment's results
        let mut msa: Vec<String> = Vec::new();
        for i in 0..n_seqs {
            let mut curr_aln = String::with_capacity(msa_l as usize);
            let outer_pointer = *msa_seq.add(i as usize);
            for j in 0..msa_l {
                let inner_pointer = *(outer_pointer.add(j as usize));
                curr_aln.push(
                    *AbpoaAligner::NT256_TABLE
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
        );

        // Read the consensus
        let mut cons = String::with_capacity(*cons_l as usize);
        for i in 0..cons_n {
            let offset = *cons_l.add((i) as usize);
            for j in 0..offset {
                let outer_pointer = *cons_seq.add(i as usize);
                let inner_pointer = *(outer_pointer.add(j as usize));
                cons.push(
                    *AbpoaAligner::NT256_TABLE
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
        /*
        let bseq: Vec<u8> = seq
            .chars()
            .map(|c| *(AbpoaAligner::NT4_TABLE).get(c as usize).unwrap())
            .collect();
        */

        let bseq: Vec<u8> = AbpoaAligner::convert_seq_to_bseq(seq);

        // First add the nodes to the graph
        // NOTE: in abpoa, each node has length 1 (i.e. a single nucleotide)
        let ids: Vec<i32> = bseq
            .into_iter()
            .map(|s| abpoa_add_graph_node((*self.ab).abg, s))
            .collect();

        assert_eq!(seq.len(), ids.len());
        for i in 0..ids.len() {
            let id = ids.get(i).unwrap();
            let seq = seq.chars().nth(i).unwrap();
            self.nodes_str.insert(*id, seq);
        }

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

    pub unsafe fn create_graph_from_paths(&mut self, seq: &str) {
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

        assert_eq!(seq.len(), ids.len());
        for i in 0..ids.len() {
            let id = ids.get(i).unwrap();
            let seq = seq.chars().nth(i).unwrap();
            self.nodes_str.insert(*id, seq);
        }

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
        let start_add_nodes = Instant::now();

        // Add nodes
        nodes.iter().for_each(|n| self.add_nodes_from_seq(n));
        info!(
            "Adding the nodes took: {} ms",
            start_add_nodes.elapsed().as_millis()
        );

        // Build a HashMap having as keys the abpoa_ids and as values the abstraction_ids
        let start_map_creation = Instant::now();
        // TODO: would like to use rayon but it does not like the .get()
        for abstraction_id in 0..self.nodes.len() {
            for abpoa_id in self.nodes.get(abstraction_id).unwrap() {
                self.abpoa_id_to_abstraction_id
                    .insert(*abpoa_id, abstraction_id);
            }
        }
        info!(
                "Creating the map took: {} ms",
                start_map_creation.elapsed().as_millis()
        );

        // Add edges between nodes
        edges.iter().for_each(|e| {
            self.edges.push((e.0, e.1));
            let last_of_start_node = self.nodes.get(e.0).unwrap().last().unwrap().clone();
            let first_of_end_node = self.nodes.get(e.1).unwrap().first().unwrap().clone();
            self.add_edge(last_of_start_node, first_of_end_node);
        });

        let start_head_tails = Instant::now();
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
            info!(
                "Finding heads and tails took: {} ms",
                start_head_tails.elapsed().as_millis()
            );
        }
    }

    /*
    // This adds an edge start->node for each abpoa_id in the abstraction nodes with no predecessors
    pub unsafe fn find_heads(&self) -> Vec<i32> {
        let abpoa_nodes_without_predecessor: Vec<i32> = self
            .nodes
            .clone()
            .into_iter()
            .flatten()
            .filter(|abpoa_node| self.is_head(abpoa_node))
            .collect();
        //println!("Abpoa nodes without predecessor: {:?}", abpoa_nodes_without_predecessor);
        //println!("Abpoa id to abstraction id: {:?}", self.abpoa_id_to_abstraction_id);

        let mut abstraction_nodes_without_predecessor: Vec<usize> = abpoa_nodes_without_predecessor
            .iter()
            .map(|node| self.abpoa_id_to_abstraction_id.get(node).unwrap().clone())
            .collect();

        abstraction_nodes_without_predecessor.sort();
        abstraction_nodes_without_predecessor.dedup();

        let mut heads: Vec<i32> = vec![];
        for abs_node in abstraction_nodes_without_predecessor {
            let abpoa_nodes = self.nodes.get(abs_node).unwrap();
            abpoa_nodes.iter().for_each(|node| heads.push(*node));
        }

        heads
    }

    // This adds an edge end->node for each abpoa_id in the abstraction nodes with no successor
    pub unsafe fn find_tails(&self) -> Vec<i32> {

        let abpoa_nodes_without_successor: Vec<i32> = self
            .nodes
            .clone()
            .into_iter()
            .flatten()
            .filter(|abpoa_node| self.is_tail(abpoa_node))
            .collect();

        let mut abstraction_nodes_without_successor: Vec<usize> = abpoa_nodes_without_successor
            .iter()
            .map(|node| self.abpoa_id_to_abstraction_id.get(node).unwrap().clone())
            .collect();

        abstraction_nodes_without_successor.sort();
        abstraction_nodes_without_successor.dedup();

        let mut tails: Vec<i32> = vec![];
        for abs_node in abstraction_nodes_without_successor {
            let abpoa_nodes = self.nodes.get(abs_node).unwrap();
            abpoa_nodes.iter().for_each(|node| tails.push(*node));
        }

        tails
    }
     */

    // This adds an edge start->node for ONLY the first abpoa_id in the abstraction_id without pred.
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

    // This adds an edge node->end for ONLY the last abpoa_id in the abstraction_id without succ.
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

    fn is_head(&self, node: &i32) -> bool {
        let node_appears_as_end_edge = self.edges_abpoa.iter().any(|x| x.1 == *node);
        return !node_appears_as_end_edge;
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

        let aln_ok = abpoa_align_sequence_to_graph(
            self.ab,
            self.abpt,
            bseq.as_mut_ptr(),
            seq.len() as i32,
            &mut res,
        );

        if aln_ok != 0 {
            panic!("Alignment did not work correctly!");
        }

        // Get cigars represented as integers
        let cigar_integers: Vec<u64> = (0..res.n_cigar)
            .map(|i| *res.graph_cigar.add(i as usize))
            .collect();
        //println!("Cigar integers len: {}", cigar_integers.len());

        // De-allocate cigar (once we have the integers, it is no longer needed)
        if res.n_cigar > 0 {
            free(res.graph_cigar as *mut c_void);
        }

        // Get cigars as records that make sense
        let cigar_vec_records_temp: Vec<AbpoaGraphCigar> = cigar_integers
            .into_iter()
            .filter_map(|curr_cigar| {
                let op = (curr_cigar & 0xf) as u32;
                //println!("Found op: {}", op);
                let mut node_id: Option<i32> = None;
                let mut query_id: Option<i32> = None;
                let mut op_len: Option<i32> = None;

                let op_char = match op {
                    ABPOA_CMATCH => {
                        // for MATCH/MISMATCH: node_id << 34  | query_id << 4 | op
                        node_id = Some(((curr_cigar >> 34) & 0x3fffffff) as i32);
                        query_id = Some(((curr_cigar >> 4) & 0x3fffffff) as i32);
                        Some('M')
                    }
                    ABPOA_CINS => {
                        // for INSERTION: query_id << 34 | op_len << 4   | op
                        query_id = Some(((curr_cigar >> 34) & 0x3fffffff) as i32);
                        op_len = Some(((curr_cigar >> 4) & 0x3fffffff) as i32);
                        Some('I')
                    }
                    ABPOA_CDEL => {
                        // for DELETION: node_id << 34  | op_len << 4   | op // op_len is always equal to 1
                        node_id = Some(((curr_cigar >> 34) & 0x3fffffff) as i32);
                        op_len = Some(((curr_cigar >> 4) & 0x3fffffff) as i32);
                        Some('D')
                    }
                    ABPOA_CDIFF => {
                        // for MATCH/MISMATCH: node_id << 34  | query_id << 4 | op
                        node_id = Some(((curr_cigar >> 34) & 0x3fffffff) as i32);
                        query_id = Some(((curr_cigar >> 4) & 0x3fffffff) as i32);
                        Some('X')
                    }
                    ABPOA_CSOFT_CLIP => {
                        // for CLIP query_id << 34 | op_len << 4   | op
                        query_id = Some(((curr_cigar >> 34) & 0x3fffffff) as i32);
                        op_len = Some(((curr_cigar >> 4) & 0x3fffffff) as i32);
                        Some('S')
                    }
                    ABPOA_CHARD_CLIP => {
                        // for CLIP query_id << 34 | op_len << 4   | op
                        query_id = Some(((curr_cigar >> 34) & 0x3fffffff) as i32);
                        op_len = Some(((curr_cigar >> 4) & 0x3fffffff) as i32);
                        Some('H')
                    }
                    _ => None,
                };

                match op_char {
                    Some(op_char) => {
                        let record =
                            AbpoaGraphCigar::cigar_from_params(op_char, node_id, query_id, op_len);
                        //println!("Found cigar record: {:?}", record);
                        Some(record)
                    }
                    _ => None,
                }
            })
            .collect();

        println!("TEMP is: {:?}", cigar_vec_records_temp.iter().map(|x| x.op).collect::<Vec<char>>());

        // Remove all the deletions on the suffix of the last abstraction node
        let mut cigar_vec_records: Vec<AbpoaGraphCigar> = vec![];
        let mut non_del_found = false;
        // Iterate over aln records until an op != 'D' is found, add everything after that
        for record in cigar_vec_records_temp.into_iter().rev() {
            if record.op != 'D' {
                non_del_found = true;
            }
            if non_del_found {
                cigar_vec_records.push(record);
            }
        }
        // Must reverse since I was adding from end to start
        cigar_vec_records.reverse();

        println!("Cigar vec records is: {:?}", cigar_vec_records.iter().map(|x| x.op).collect::<Vec<char>>());
        // Get the ids and convert them to abstraction ids
        let abpoa_ids: Vec<i32> = cigar_vec_records.iter().filter_map(|x| x.node_id).collect();
        //println!("Abpoa ids: {:#?}", abpoa_ids);

        let graph_ids: Vec<usize> = abpoa_ids
            .iter()
            .filter_map(|id| self.abpoa_id_to_abstraction_id.get(id))
            .map(|id| *id)
            .collect();
        //println!("Graph ids: {:#?}", graph_ids);

        let mut first_node_offset: usize = 0;
        let mut last_node_offset: usize = 0;
        if res.n_cigar > 0 {
            // Find start and end offset (relative to abstraction nodes)
            let first_abstraction_node: usize = *graph_ids.first().unwrap();
            let first_abpoa_node: i32 = *abpoa_ids.first().unwrap();
            let first_node_ids: Vec<i32> = self.nodes.get(first_abstraction_node).unwrap().clone();
            first_node_offset = first_node_ids.iter().position(|id| *id == first_abpoa_node).unwrap();

            let last_abstraction_node: usize = *graph_ids.last().unwrap();
            let last_abpoa_node: i32 = *abpoa_ids.last().unwrap();
            let last_node_ids: Vec<i32> = self.nodes.get(last_abstraction_node).unwrap().clone();
            last_node_offset = last_node_ids.iter().position(|id| *id == last_abpoa_node).unwrap();
        }


        /*
        let new_cigar_string : Vec<Option<AbpoaGraphCigar>> =
            once(None)
                .chain(cigar_vec_records.into_iter().map(|x| Some(x)))
                .chain(once(None))
                .into_iter().windows(2);
         */

        //let cigar_vec: Vec<char> = cigar_vec_records.iter().map(|rec| rec.op).collect();
        let mut cigar_vec: Vec<char> = Vec::new();

        // Get the cigar as a string
        let mut new_cigar_string = String::new();
        let mut last_cigar_val: Option<AbpoaGraphCigar> = None;
        let mut last_cigar_count: u64 = 0;
        for (i, curr_cigar) in cigar_vec_records.iter().enumerate() {
            match last_cigar_val {
                // First iteration
                None => (),
                // Second to last iteration
                Some(last_val) => match (curr_cigar.op != last_val.op, last_val.op) {
                    (true, 'I') => {
                        let ins_len = last_val.op_len.unwrap();
                        new_cigar_string.push_str(&mut format!("{}{}", ins_len, last_val.op));

                        for i in 0..ins_len {
                            cigar_vec.push(last_val.op);
                        }

                        last_cigar_count = 0;
                    }
                    (true, _) => {
                        new_cigar_string
                            .push_str(&mut format!("{}{}", last_cigar_count, last_val.op));

                        for i in 0..last_cigar_count {
                            cigar_vec.push(last_val.op);
                        }

                        last_cigar_count = 0;
                    }
                    (false, _) => (),
                },
            }

            last_cigar_val = Some(curr_cigar.clone());
            last_cigar_count += 1;

            // Last value check
            if i == cigar_vec_records.len() - 1 {
                let last_val = last_cigar_val.unwrap();
                match curr_cigar.op {
                    'I' => {
                        let ins_len = last_val.op_len.unwrap();
                        new_cigar_string.push_str(&mut format!("{}{}", ins_len, last_val.op));

                        for i in 0..ins_len {
                            cigar_vec.push(last_val.op);
                        }

                        last_cigar_count = 0;
                    }
                    _ => {
                        new_cigar_string
                            .push_str(&mut format!("{}{}", last_cigar_count, last_val.op));

                        for i in 0..last_cigar_count {
                            cigar_vec.push(last_val.op);
                        }

                        last_cigar_count = 0;
                    }
                }
            }
        }

        // Obtain cs string
        let mut cs_string = String::new();
        /*
        let mut match_count = 1;
        let mut tmp_string: String = String::new();
        let mut last_char = ' ';
        let mut char = ' ';

        for i in 0..cigar_vec.len() {
            char = *cigar_vec.get(i).unwrap();

            if char != last_char {
                if last_char != ' ' {
                    // This is the initial delimiter for this loop
                    match last_char {
                        'M' => cs_string.push_str(&mut format!("{}", match_count)),
                        'I' => cs_string.push_str(&mut format!("{}{}", '-', tmp_string)),
                        'D' => cs_string.push_str(&mut format!("{}{}", '+', tmp_string)),
                        _ => (),
                    }
                }
                last_char = char;
                // Reset tmp variables
                match_count = 1;
                tmp_string = String::new();
            } else {
                match char {
                    'M' => match_count += 1,
                    // TODO: fix out of bounds here
                    'I' => tmp_string.push(seq.char_indices().nth(i).unwrap().1),
                    'D' => {
                        let id = abpoa_ids.get(i).unwrap();
                        tmp_string.push(*self.nodes_str.get(id).unwrap())
                    }
                    _ => (),
                }
            }
        }

        // Last iteration
        match char {
            'M' => cs_string.push_str(&mut format!("{}", match_count)),
            'I' => cs_string.push_str(&mut format!("{}{}", '-', tmp_string)),
            'D' => cs_string.push_str(&mut format!("{}{}", '+', tmp_string)),
            _ => (),
        }

        //assert_eq!(abpoa_ids.len(), graph_ids.len());
         */

        AbpoaAlignmentResult::new_with_params(
            new_cigar_string.as_str(),
            abpoa_ids,
            graph_ids,
            first_node_offset,
            last_node_offset,
            &res,
            cigar_vec,
            cs_string.as_str(),
        )
    }

    pub unsafe fn create_align_safe(
        nodes: &Vec<&str>,
        edges: &Vec<(usize, usize)>,
        query: &str,
        mode: AbpoaAlignmentMode,
    ) -> AbpoaAlignmentResult {
        // Create the aligner
        let mut aligner = AbpoaAligner::new_with_params(6,9,10,mode);

        // Create the graph
        aligner.add_nodes_edges(nodes, edges);

        // Align sequence to graph
        let res = aligner.align_sequence(query);

        // IMPORTANT: delete memory associated with ab and abpt
        abpoa_free(aligner.ab);
        abpoa_free_para(aligner.abpt);

        res
    }

    pub unsafe fn read_graph_from_file(&mut self, path_file: &PathBuf) {
        let parser = GFAParser::new();
        let gfa: GFA<usize, ()> = parser.parse_file(&PathBuf::from(path_file)).unwrap();
        self.read_graph_from_GFA(&gfa);
    }

    pub unsafe fn read_graph_from_GFA(&mut self, gfa: &GFA<usize, ()>) {
        let nodes_as_strings: Vec<String> = gfa
            .segments
            .iter()
            .map(|s| s.sequence.to_string())
            .collect();
        let nodes: Vec<&str> = nodes_as_strings.iter().map(|s| s.as_str()).collect();
        let edges: Vec<(usize, usize)> = gfa
            .links
            .iter()
            .map(|e| (e.from_segment - 1, e.to_segment - 1))
            .collect();

        //println!("Nodes from GFA: {:?}", nodes_as_strings);
        //println!("Edges from GFA: {:?}", edges);

        self.add_nodes_edges(&nodes, &edges);
    }
}

#[derive(Debug, Clone, Copy)]
pub struct AbpoaGraphCigar {
    op: char,
    node_id: Option<i32>,
    query_id: Option<i32>,
    op_len: Option<i32>,
}

impl AbpoaGraphCigar {
    pub fn cigar_from_params(
        op: char,
        node_id: Option<i32>,
        query_id: Option<i32>,
        op_len: Option<i32>,
    ) -> Self {
        AbpoaGraphCigar {
            op,
            node_id,
            query_id,
            op_len,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::abpoa::{abpoa_generate_gfa, stdout};
    use gfa::gfa::Line::Path;

    #[test]
    fn test_simply_create() {
        // to be run with: (requires rust nightly -- rustup default nighly)
        // RUSTFLAGS="-Z sanitizer=address" cargo test test_simply_create --target x86_64-unknown-linux-gnu
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&vec!["AAA", "CC"], &vec![(0, 1)]);
            aligner.align_sequence("AACC");
            // if these two functions are not called, there's a leak
            abpoa_free(aligner.ab);
            abpoa_free_para(aligner.abpt);
        }
    }

    #[test]
    fn test_create_and_align() {
        // to be run with: (requires rust nightly -- rustup default nighly)
        // RUSTFLAGS="-Z sanitizer=address" cargo test test_create_and_align --target x86_64-unknown-linux-gnu
        unsafe {
            let result = AbpoaAligner::create_align_safe(&vec!["AAA", "CC"], &vec![(0, 1)], "AACC",AbpoaAlignmentMode::Global);
        }
    }

    /*
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
     */

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
            println!("Res: {:?}", res);

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
            println!("Res: {:?}", res);

            //Unexpected but same behavior in C ver
            assert_eq!(res.cigar, String::from("3M"));
            // These are the nodes in abpoa (remember, nodes are 1-base only!)
            assert_eq!(res.abpoa_nodes, vec![2, 3, 4]);
            // These are the nodes in our graph abstraction (nodes can have length > 1)
            assert_eq!(res.graph_nodes, vec![0, 0, 0]);
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
            let mut abpoa_ids: Vec<i32> = Vec::new();
            let mut graph_ids: Vec<usize> = Vec::new();
            let mut cigar_vec: Vec<char> = Vec::new();

            // Navigate the cigar
            let mut op: u32 = 0;
            let mut op_char = ' ';
            let mut node_id: i32 = 0;
            for i in 0..res.n_cigar {
                let curr_cigar = res.graph_cigar.add(i as usize);
                op = (*curr_cigar & 0xf) as u32;

                // TODO: controlla S e H
                // usa abpoa da linea di comando

                op_char = match op {
                    ABPOA_CMATCH => 'M',
                    ABPOA_CINS => 'I',
                    ABPOA_CDEL => 'D',
                    ABPOA_CDIFF => 'X',
                    ABPOA_CSOFT_CLIP => 'S',
                    ABPOA_CHARD_CLIP => 'H',
                    _ => ' ',
                };

                node_id = ((*curr_cigar >> 34) & 0x3fffffff) as i32;
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
                0,
                0,
                &res,
                cigar_vec,
                String::new().as_str(),
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
            let mut abpoa_ids: Vec<i32> = Vec::new();
            let mut graph_ids: Vec<usize> = Vec::new();
            let mut cigar_vec: Vec<char> = Vec::new();

            // Navigate the cigar
            let mut op: u32 = 0;
            let mut op_char = ' ';
            let mut node_id: i32 = 0;
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

                node_id = ((*curr_cigar >> 34) & 0x3fffffff) as i32;

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
                0,
                0,
                &res,
                cigar_vec,
                String::new().as_str(),
            )
        }
    }

    #[test]
    fn add_complex_graph_1() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&vec!["ACG", "AAA"], &vec![(0, 1)]);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
        }
    }

    #[test]
    fn add_complex_graph_2() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&vec!["ACG", "AAA", "CC"], &vec![(0, 1)]);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
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
            println!("Result is: {:#?}", result);

            assert_eq!(
                result.abpoa_nodes,
                vec![
                    2, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                    27, 28, 29, 30
                ]
            );

            assert_eq!(
                result.graph_nodes,
                vec![0, 2, 2, 2, 2, 4, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 9, 9, 9, 9]
            );
        }
    }

    #[test]
    fn safe_function_test() {
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

            let res = AbpoaAligner::create_align_safe(&nodes, &edges, "AGAAT", AbpoaAlignmentMode::Global);
            println!("{:#?}", res.cigar)
        }
    }

    #[test]
    fn test_cigar_length() {
        unsafe {
            let nodes = vec!["ACGT", "TTG", "CGA"];
            let edges = vec![(0, 1), (1, 2)];
            let query = "ACTTTGCGTTTTTTT";

            let total_node_length: usize = nodes.iter().map(|x| x.len()).sum();

            let res = AbpoaAligner::create_align_safe(&nodes, &edges, &query, AbpoaAlignmentMode::Global);
            println!("{:#?}", res.cigar_vec);
            //assert_eq!(total_node_length, res.cigar_vec.len())
        }
    }

    #[test]
    fn test_cigar_length_suffix() {
        unsafe {
            let nodes = vec!["ACGT", "TTG", "CGA"];
            let edges = vec![(0, 1), (1, 2)];
            let query = "TGCG";

            let total_node_length: usize = nodes.iter().map(|x| x.len()).sum();

            let res = AbpoaAligner::create_align_safe(&nodes, &edges, &query, AbpoaAlignmentMode::Global);
            println!("{:#?}", res.cigar_vec);
            //assert_eq!(res.cigar_vec.len(), 5)
        }
    }

    #[test]
    fn test_validate_with_GFA() {
        unsafe {
            let nodes = vec!["ACGT", "TTG", "CGA"];
            let edges = vec![(0, 1), (1, 2)];

            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&nodes, &edges);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
            println!("");

            let query = "ACTTTGCGTTTTTTT";
            let result = aligner.align_sequence(query);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
        }
    }

    #[test]
    fn test_validate_with_GFA_2() {
        unsafe {
            let nodes = vec!["TTG", "G", "AAT"];
            let edges = vec![(0, 1), (1, 2)];

            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&nodes, &edges);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
            println!("");

            let query = "AAT";
            let result = aligner.align_sequence(query);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
            println!("Result: {:#?}", result);
        }
    }

    #[test]
    fn test_validate_0() {
        // Funziona come abpoa (grafo con path)
        unsafe {
            let nodes = vec!["ACG"];
            let edges = vec![];

            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&nodes, &edges);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
            //println!("");
            //println!("");

            let query = "GTT";
            let result = aligner.align_sequence(query);

            //abpoa_Fgenerate_gfa(aligner.ab, aligner.abpt, stdout);
            //println!("Result: {:#?}", result);
        }
    }

    #[test]
    fn test_validate_1() {
        // Funziona come abpoa (grafo con path)
        unsafe {
            let nodes = vec!["ACGA", "TTG", "CGA"];
            let edges = vec![(0, 1), (1, 2)];

            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&nodes, &edges);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
            //println!("");

            let query = "ACGA";
            let result = aligner.align_sequence(query);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
            println!("Result: {:#?}", result);
        }
    }

    #[test]
    fn test_validate_2() {
        unsafe {
            let nodes = vec!["ACGT", "TTG", "AAA", "CGA"];
            let edges = vec![(0, 1), (0, 2), (1, 3), (2, 3)];

            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&nodes, &edges);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
            //println!("");

            let query = "ACGTAAACGATTT";
            let result = aligner.align_sequence(query);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
            println!("Result: {:#?}", result);
        }
    }

    #[test]
    fn test_validate_3() {
        unsafe {
            let nodes = vec!["ACGT", "TTG", "AAA", "CGA"];
            let edges = vec![(0, 1), (0, 2), (1, 3), (2, 3)];

            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&nodes, &edges);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
            //println!("");

            let query = "TTTTTTG";
            let result = aligner.align_sequence(query);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
            println!("Result: {:#?}", result);
        }
    }

    #[test]
    fn test_validate_subgraph_7() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();
            let nodes = vec![
                "TTG",
                "A",
                "G",
                "AAAT",
                "AA",
                "TTTCT",
                "GG",
                "AGTTCTAT",
                "A",
                "T",
                "ATAT",
                "A",
                "T",
                "CCAACTCTCTG",
            ];
            let edges = vec![
                (0, 1),
                (0, 2),
                (1, 3),
                (2, 3),
                (3, 4),
                (3, 5),
                (4, 5),
                (5, 6),
                (5, 7),
                (6, 7),
                (7, 8),
                (7, 9),
                (8, 10),
                (9, 10),
                (10, 11),
                (10, 12),
                (11, 13),
                (12, 13),
            ];
            aligner.add_nodes_edges(&nodes, &edges);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
            //println!("Edges: {:#?}", aligner.edges_abpoa);
            //println!("Nodes: {:#?}", aligner.abpoa_id_to_abstraction_id);

            /*
            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.add_nodes_edges(&nodes, &edges);
            abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
            //println!("");
            */

            let query = "CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG";
            let result = aligner.align_sequence(query);
            println!("Result: {:#?}", result);
        }
    }

    #[test]
    fn test_validate_subgraph_8() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();
            let nodes = vec![
                "TTG",
                "A",
                "G",
                "AAAT",
                "AA",
                "TTTCT",
                "GG",
                "AGTTCTAT",
                "A",
                "T",
                "ATAT",
                "A",
                "T",
                "CCAACTCTCTG",
            ];
            let edges = vec![
                (0, 1),
                (0, 2),
                (1, 3),
                (2, 3),
                (3, 4),
                (3, 5),
                (4, 5),
                (5, 6),
                (5, 7),
                (6, 7),
                (7, 8),
                (7, 9),
                (8, 10),
                (9, 10),
                (10, 11),
                (10, 12),
                (11, 13),
                (12, 13),
            ];
            aligner.add_nodes_edges(&nodes, &edges);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
            //println!("Edges: {:#?}", aligner.edges_abpoa);
            //println!("Nodes: {:#?}", aligner.abpoa_id_to_abstraction_id);

            let query = "CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG";
            let result = aligner.align_sequence(query);
            //abpoa_generate_gfa(aligner.ab, aligner.abpt, stdout);
            println!("Result: {:#?}", result);
        }
    }

    #[test]
    fn test_read_graph() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.read_graph_from_file(&PathBuf::from("test/graph.gfa"));
            //println!("Nodes: {:#?}", aligner.nodes);
            //println!("Edges: {:#?}", aligner.edges);
        }
    }

    #[test]
    fn test_with_luca() {
        unsafe {
            let mut aligner = AbpoaAligner::new_with_example_params();
            aligner.read_graph_from_file(&PathBuf::from("test/graph.gfa"));
        }
    }

    #[test]
    fn test_offsets_aln() {
        unsafe {
            let nodes = vec!["ACGT", "TTG", "CGA"];
            let edges = vec![(0, 1), (1, 2)];
            let query = "CGTTTGCG";

            let total_node_length: usize = nodes.iter().map(|x| x.len()).sum();

            let res = AbpoaAligner::create_align_safe(&nodes, &edges, &query, AbpoaAlignmentMode::Global);
            assert_eq!(res.aln_start_offset, 1);
            // TODO: fix this
            //assert_eq!(res.aln_end_offset, 1);
        }
    }

    #[test]
    fn test_offsets_aln_2() {
        unsafe {
            let nodes = vec!["ACGTTTTCGA"];
            let edges = vec![];
            let query = "TTTTCGA";

            let total_node_length: usize = nodes.iter().map(|x| x.len()).sum();

            let res = AbpoaAligner::create_align_safe(&nodes, &edges, &query, AbpoaAlignmentMode::Global);
            assert_eq!(res.aln_start_offset, 3);
            assert_eq!(res.aln_end_offset, 9);
        }
    }

    #[test]
    fn figure_out_whats_wrong() {
        unsafe {
            let nodes = vec!["CTCCTTCTTGGGCTAGGACTGTGCCCACAGCTGACAGACCTCAAACAGTAGAAGAAACAGGGATGGAGGCCAGAATACCACTCCTCCCTTGGATCAGGAGAGGGAGCTGTCACCTGAGGTACAGGAGATCCTATACCACAGAGTGACTCTCTTAAAGGGCCAGACCTCTCTCAGGGGCAATTAAGGAATCTAGTCTCGCTGGAGATTCCATCCTTCAGATGAACTGATGAGCAGTTCTCTTTGACTCCCAGTATTAGGAATCACGGGGGAGTTTCTCTCGTGCCTGATTCTCAGCCCCACACCAAGAGTTTTTGGAGGTCTGACTCCAGCTTTTCTCAGTCACTCAGCATCCACACAGGCCAGGACCAGAAATCCCTTTTCACCTTCTACCCTGGGCTAGCTCATCCCGATTCTAGAACTTTCCAAGGAATAAGAGGCTATCCCAGATCCCTAAGTCCAGGCTGGTGTCAAGGTTTTGTCCTCTTCTCCTACTATAATTGTCCTCTTCCTTCTCAGGATGGTCACATGGGTGCTGCTGGAGTGTCCCATGAGAGATACAAAGTGCCTGAATTTTCTGACTCTTCCCCTCAGAGCCCCCAAAGACACACGTGACTCACCACCCCATCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCCCTGGGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCAGGATGGGGAGGGCCATACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCAGGGGATGGAACCTTCCAGAAGTGGGCAGCTGTGGTGGTGCCTTCTGGAGAGGAGCAGAGATACACGTGCCATGTGCAGCATGAGGGGCTACCCGAGCCCGTCACCCTGAGATGGAGTAAGGAGGGGGATGGGAGGTCATGTCTCTTCTCAGGGAAAGCGGGAGCCCTTCTGGAGCCCTTCCGCAGGGTCAGGGCTGAGGCCTGGGGGTCAGGGCCCCTTACGTTCCCCTCTTTTCCCAGAGCCGGCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCCTTGGATCTGTGGTCTCTGGAGCTGTGGTTGCTGCTGTGATATGGAGGAAGAAGAGCTCAGGTGGGGAAGGGAGAAGGGTGGGGTCTGAGTTTTCTTGTCCCACTGGGTGTTTCAAGCCCTAGGTAAAAGTGTGTCCTGCCTCGTTACTGGGAAGCACCATCCACACACACGAGCCTACCCAGCCTGGGGCCCTGTGTGCCAGCACCTACTCTTTTTTTTTGAGACGGAGTCTTGGCTCTGTCACCCAGGCTGGAGTGCAATGGCGTGGTTTCAGCTCACTGCAACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTCAGCCTCCCTAGTAGCTGGGACTACACATGCGTGCCACCACACCTGGCTAATTTTTTTTTTTGTATTTTTAGTGGAGATGGGGTTTCACTATGTTGGCCAGGCTGGTCTCGAACTCCTGACTTTGTGATCTGCCTGCCTCGGCCTCCCAAAGTGCTGGGATTACAGTCGTGAGCCACCGCACCCAGCCGCACCTACTCTTTTGTAAAGCACCTGTGACAATGAAGGACAGATTTATCACCTTGACGATTGTGGTGATGGGGACCTGATCCCAGCAGTCACAGGTCACAGGGGAAGGTCCCTGCTGAAGACAGACCTCAGAAGGGCAGTTGATCCAGGACCCACACCTGCTTTCTTCACGTTTCCTGATCCTGCCCTGGGTCTGCAGTCACAGTTCAGGAAACTTCTCTGGGATCCAAAACTAGGAGGTTCCTCTAGGACCTTATGGCCCTGCCTCCTCCCTGGCCCCTCACAGGACATTTTCTTCCAACAGGTGGAAAAGGAGGGAGCTACTCTAAGGCTGAGTGTAAGTGCGGGGCGGGAGCGTGGAGGAGCTCGCCCACCCTATAATTCCTCCTGCACCACATCTCCTGTGGGCTCTGACCAGGTCTTGTTTTTGTTCTACCCCAGGGAGCGACAGTGCCCAGGGGTCTGAGTCTCACAGCTTGTAAAGGTGAGATTCTGGGGGTCTGAAGTGGGTGGAGGGTGGGGCAGAGGGGACAGGACTGGGTTGTGGGGATTTTTTGATTCAGAATTTTTGAGTGTGTGGTGGGCTGTTCAGAGTGTCATCACTTACCGTGACTGACCTGAATTTGTTCATGACTATTTTCTTCTGTAGCCTGAGACAGCTGCCTTGTGTGCGACTGAGATGCACAGCTGCCTTGTGTGCGACTGAGATGCAGGATTTCCTCACGCCTCCCCTATGTGTCTTAGGGGACTCTGGCTTCTCTTTTTGCAAGGGCCTCTGAATCTGTCTGTGTCCCTGTTAGCACAATGTGAGGAGGTAGAGAAACAGTCCACCTCTGTGTCTACCATGACCCCCTTCCTCACACTGACCTGTGTTCCTTCCCTGTTCTCTTTTCTATTAAAAATAAGAACCTGGGCAGAGTGCGGCAGCTCATGCCTGTAATCCCAGCACTTAGGGAGGCCGAGGAGGGCAGATCACGAGGTCAGGAGATCGAAACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAAATACAAAAAATTAGCTGGGCGCAGAGGCACGGGCCTGTAGTCCCAGCTACTCAGGAGGCGGAGGCAGGAGAATGGCGTCAACCCGGGAGGCGGAGGTTGCAGTGAGCCAGGATTGTGCGACTGCACTCCAGCCTGGGTGACAGGGTGAAACGCCATCTCAAAAAATAAAAATT"];
            let edges = vec![];
            let query = "AAAAAATACAAAAAATTAGCTGGGCGCAGAGGCACGGGCCTGTAGTCCCAGCTACTCAGGAGGCGGAGGCAGGAGAATGGCGTCAACCCGGGAGGCGGAG";

            let res = AbpoaAligner::create_align_safe(&nodes, &edges, &query, AbpoaAlignmentMode::Global);
            println!("res is: {:?}", res);
        }
    }
}
