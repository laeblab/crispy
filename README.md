* Description:
  * crispy++ has 2 functions
  1.  count:	(indexing) find and count all PAM adjacent sequences in genome.
  2.  score:	Use index file (counts.laeb) to score list of supplied targets

* Enzyme support:
  1.  spCas9
  2.  (MAD7, not ready yet)

* Quick start:
  1. Compile by typing "make" in root directory
  2. Run example by typing this in root directory
```bash
./bin/crispy count Cas9 test_data/test_genome_fasta.fa
diff cas9_counts.laeb test_data/cas9_counts.laeb
./bin/crispy score 8 test_data/test_targets.tsv cas9_counts.laeb
diff scored.laeb test_data/scored.laeb
```
If the diff commands returned nothing, then everything works.
You can now run "count" with your genome of choice and then construct a target file to run "score" on.

* Target file format:
  * A tab separated file.
  * No header.
  * Each line is a target. Only first column is used, everything after is left as is in output file.
  * First column must be a 20 bp DNA sequence. E.g. you want to score a Cas9 target ACGTACGACGAATCGATCGA-NGG, then just write in the first 20 bp, that is the whole sequence, up to, but not including, the PAM sequence.
  * See test_data/test_targets.tsv for a valid target file

* Output file format [scored.laeb]
  * Identical to targets tsv file, except 1 column has been added which is the score (lower is better, but be aware that invalid lines will also get a 0)
  
* Dev instructions:
  * test_data folder contains a test_fasta.fa and 2 .laeb files. If you make changes to count.cpp, please make sure that the output is identical (diff should give no output) with the relevant .laeb file.