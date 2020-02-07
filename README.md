# CRISPyR - A CRISPR target finding tool based on CRISPy

CRISPyR is a command-line tool based on CRISPy ([https://doi.org/10.1002/bit.25233](https://doi.org/10.1002/bit.25233)).
Its a RUST rewrite (and extension) of the original python code, hence the 'R' in the name.
You can use it for finding and scoring CRISPR target sites in FASTA sequences. 
CRISPyR supports the Cas9 (the default) and the Mad7 (experimental) endonucleases.

CRISPyR finds candidate target sequences by searching a genome for PAM sites
(NGG-'3 for Cas9, 5'-YTTN for Mad7). Alternatively, the user may supply a list
of target sequences for a given enzyme. Each target sequence is scored for off-
targets based on the number of sequences matching the 13bp adjacent to the PAM
(upstream of the PAM for Cas9 and downstream for Mad7). The resulting score
represents the number putative off-targets, with a higher score representing
more and/or closely matching off-target sequences.

CRISPyR can return also return a list of all positions in the genome that it
considers as putative off-targets for a given target sequence.

If you use CRISPyR, then please cite the paper:

    Ronda et al (2014). Accelerating genome editing in CHO cells using CRISPR
    Cas9 and CRISPy, a web-based target finding tool.Biotechnol. Bioeng. 111:
    1604–1616. doi: 10.1002/bit.25233

[https://doi.org/10.1002/bit.25233](https://doi.org/10.1002/bit.25233)


## Installation

CRISPyR is built using the [Rust](https://www.rust-lang.org) programming
language. To compile CRISPyR, install Rust as described on the Rust website
and run 'cargo build --release' in the source directory:

    $ git clone https://github.com/laeblab/crispy.git
    $ cd crispy
    $ cargo build --release

The resulting executable is located at './target/release/crispyr'.


## Usage

The following examples use the files included in the 'examples' folder.

Note that CRISPyR outputs results as tab separated tables. The results shown
below have been piped through 'column -t' for the sake of readability.


### Indexing a genome

To use CRISPyR, you must first index the genome against which you wish to
score potential guide RNAs. To index a genome, use the CRISPyR 'index' command
on a FASTA file:

    $ crispyr index examples/genome.fasta

This will create the index file 'examples/genome.fasta.crispyr_cas9'.

By default CRISPyR will search for PAM sites for CAS9 (NGG), but MAD7 (YTTN)
is also supported, and may be selected using the --enzyme option:

    $ crispyr index --enzyme mad7 examples/genome.fasta

This will create the index file 'examples/genome.fasta.crispyr_mad7'. Note that
the cut-site for Mad7 is not documented, and results will therefore not show
the expected cut-site for target sequences.

The index command can also be run with '--positions' to record the positions of
all PAM sites in the genome. This takes longer and significantly increases the
size of the index, but is required to run the 'offtargets' command.


### Finding target sequences

The 'find' command takes a CRISPyR index and a FASTA file as input, and prints
a table of target sites in the FASTA file along with the CRISPy score
representing the matching 13bp kmers found in the index (see above):

    $ crispyr find examples/genome.fasta.crispyr_cas9 examples/genome.fasta
    Sequence                 Name   Start  End  Cutsite  Strand  Score
    GCTAGCTAGCTAGCTATAAAagg  test   14     36   31       +       1000
    CTAGCTAGCTAGCTATAAAAggg  test   15     37   32       +       1000
    CTACGTAGCTACTAGCTGACtgg  test   49     71   66       +       1000
    [...]

Positions are given using base-1 values on the forward strand (regardless of
which strand the target was found). The PAM is indicated in the sequence using
lower-case letters.

A higher Score indicates more/better matching off-targets, meaning that target
sequences with lower scores should be selected when possible.


### Scoring existing target sequences

The 'score' command can be used to score a set of gRNAs that have been obtained 
without using CRISPyR. It needs an index and a tab separated table and it will
append the CRISPy score for each row in which the first column contains a
target sequence matching the PAM used to create the index:

    $ crispyr score examples/genome.fasta.crispyr_cas9 examples/targets.tsv
    AATAGCTAGCTAGCTATAAAAGG  2 copy target                      1000
    CAGCTACTAGCTAGTCGATGNGG  2 copy target                      1000
    CCCCCCCCCCCCCCCCCCCCCAA  not a target at all. should get 0  0

CRISPyR will only attempt to score target sequences that 
1. contains a PAM,
2. contains at least 13 bp in addition to the PAM
3. the 13 bp only consist of nucleotides A, C, G, or T. Invalid target sequences will be assigned the
score 'NA'.

The input table may (optional) contain a header, in which case CRISPyR will
appends a column named 'Score'.


### Finding potential off-targets

The 'offtargets' command takes a CIRPSyR index and a tab separated table and
returns a list of sites for each (unique) target sequence that CRISPyR considers
to be off-targets for that sequence:

    $ crispyr offtargets examples/genome.fasta.crispyr_cas9 examples/targets.tsv
    Query                    Offtarget                Name   Start  End  Cutsite  Strand  Score
    aatagctAGCTAGCTATAAAagg  gctagctAGCTAGCTATAAAagg  test   14     36   31       +       500
    aatagctAGCTAGCTATAAAagg  gctagctAGCTAGCTATAAAagg  test1  14     36   31       +       500
    cagctacTAGCTAGTCGATGngg  cagctacTAGCTAGTCGATGcgg  tesss  120    142  137      +       500
    cagctacTAGCTAGTCGATGngg  cagctacTAGCTAGTCGATGcgg  tesss  162    184  179      +       500

The 13bp kmer used to identify off-targets (see above) is written in uppercase.

The Score column represents the weight given to each individual off-target, with
the score of the query being the sum of those scores (see src/score.rs). It is
also possible to filter off-targets based on this score, with a score of 500
representing perfect matches between the 13bp kmer in the query and off-target:

    $ crispyr offtargets examples/genome.fasta.crispyr_cas9 examples/targets.tsv --min-score 500

For the 'offtargets' command to work, the genome must have been indexed with the
'--positions' option. Note that this greatly increases the size of the index and
isn't recommended unless you specifically need the 'offtargets' functionality:

    $ crispyr index --positions examples/genome.fasta

Target sequences must meet the requirements described for the 'score' command;
any value that does not meet these requirements trigger a warning.

The Offtarget column requires the availability of a faidx indexed FASTA file:

    $ samtools faidx /path/to/genome.fasta

By default CRISPyR will attempt to the filename produced by removing the
'.crispyr_*' extension from the index filename, but an alternative path may be
specified using the '--fasta' option. If no such FASTA file is available, then
this column will contain the value 'NA'.
