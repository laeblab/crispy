use crate::common::KMer;
use crate::constants::*;
use crate::index::{KMerIndex, Position};

pub const MAX_MUTATIONS: usize = 3;
pub const MAX_SEED_MUTATIONS: usize = 2;
pub const SEED_SIZE: usize = 5;

// Mismatch offtarget scoring matrix
// seed mismatches on column and non-seed mismatches on row:
// 		0						1						2					3				4
// 0	0 seed and 0 non-seed	0 seed and 1 non-seed	0 seed 2 non-seed	0 seed 3 non    0 seed 4 non
// 1	1 seed and 0 non-seed	1 seed and 1 non-seed	1 seed 2 non-seed	1 seed 3 non	not-a-target
// 2	2 seed and 0 non-seed	2 seed and 1 non-seed	2 seed 2 non-seed	not-a-target	not-a-target

// A gRNA is scored by summing the score of all possible off-targets, using the above matrix.
// For example, a gRNA might have just 1 potential off-target, with 0 mismatches in the seed
// region and 0 mismatches in the rest of the sequence.assert. This gRNA would receive a score
// of 500 (1 * matrix[0][0]). Another gRNA might have 499 possible off-targets, each with 2
// mismatches in the seed region and 2 mismatches in the rest of the sequences. That gRNA would
// receive a score of 499 (499 * matrix[2][2]), making it a better choice than the first gRNA.

// Depending on your use-case you might want to change this. E.g. you could imagine you wanted
// to perform a test of gRNA edits at all offtargets and for that reason you wanted as few
// potential offtargets as possible. Then you might set the score of all situations to 1 and
// then just pick the gRNA with the lowest score.
static SCORE_MATRIX: [[u64; 5]; 3] = [[500, 100, 50, 20, 3], [80, 30, 15, 2, 0], [20, 5, 1, 0, 0]];

struct Permutation {
    kmer: KMer,
    n_seed: usize,
    n_rest: usize,
}

impl Permutation {
    fn score(&self) -> u64 {
        SCORE_MATRIX[self.n_seed][self.n_rest]
    }
}

fn permute(
    permutations: &mut Vec<Permutation>,
    kmer: KMer,
    seed_muts: usize,
    rest_muts: usize,
    pos: usize,
) {
    if seed_muts + rest_muts < MAX_MUTATIONS {
        for pos in pos..KMER_LEN {
            let (seed_muts, rest_muts) = if pos < SEED_SIZE {
                (seed_muts + 1, rest_muts)
            } else {
                (seed_muts, rest_muts + 1)
            };

            if seed_muts <= MAX_SEED_MUTATIONS {
                let current_nucleotide = (kmer.0 >> (2 * pos)) & 3;
                let kmer_template = kmer.0 & !(kmer.0 & (3 << (2 * pos)));
                for nucleotide in 0..4 {
                    if nucleotide != current_nucleotide {
                        let new_kmer = KMer::new(kmer_template | (nucleotide << (2 * pos)));
                        permutations.push(Permutation {
                            kmer: new_kmer,
                            n_seed: seed_muts,
                            n_rest: rest_muts,
                        });

                        permute(permutations, new_kmer, seed_muts, rest_muts, pos + 1);
                    }
                }
            }
        }
    }
}

fn permutations(kmer: KMer) -> Vec<Permutation> {
    let mut permutations = Vec::with_capacity(9 * 1024);
    permutations.push(Permutation {
        kmer,
        n_seed: 0,
        n_rest: 0,
    });

    permute(&mut permutations, kmer, 0, 0, 0);
    permutations
}

pub fn calculate_score(index: &KMerIndex, kmer: KMer) -> u64 {
    let mut score = 0;
    for permutation in permutations(kmer) {
        if let Some(count) = index.get_count(permutation.kmer) {
            score += u64::from(count) * permutation.score();
        }
    }

    score
}

pub fn find_offtargets(index: &KMerIndex, kmer: KMer, min_score: u64) -> Vec<(u64, &Position)> {
    let mut result = Vec::new();

    for permutation in permutations(kmer) {
        let score = permutation.score();

        if score >= min_score {
            if let Some(positions) = index.get_positions(permutation.kmer) {
                for pos in positions {
                    result.push((score, pos));
                }
            }
        }
    }

    result
}
