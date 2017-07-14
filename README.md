# PolyFud-BPP

This repository contains a prototype of PolyFud-BPP which uses cis-elements of polyadenylation sites(PolyFud) and 3'UTR base pair probabilities to predict functional polyadenylation motifs.

Requirements to compile the project: 
- Boost.Filesystem, Boost.Spirit, Boost.Foreach, Boost.Optional (used was version 1.58, Boost.Filesystem needs a compiled library)
- SeqAn's FASTA index (used was version 2.0.0)
- ViennaRNA (used was version 2.23)

The program also needs the human reference genome. There is a script in the scripts folder that will download it automatically and place it in "bin/reference\_genome/hg19" (requires max 7GB disk space). The chromosome FASTAs can be deleted afterwards.

Folder content:
- The folder perf\_testing contains the source code for the binaries that compute the sensitivity and specificity. It also contains the necessary dataset positiveSet.fa (PASes with their +/-250nt long flanking regions), negativeSet.fa (randomized sequences of positiveSet.fa), utrBppPerTranscript.txt (BPPs of positive dataset), utrBppPerTranscriptTn.txt (BPPs of negative dataset). 
- The src folder contains the source code for the refGene parser (reGene.txt in bin/ucsc\_data folder), PolyFud and PolyFudBpp. 
- unit\_tests contains the unit tests for PolyFud, PolyFudBpp and the polar utility functions.

Notices:
- Parameters are hard coded and need to be changed inside the source files! (TODO: a config file)
- PolyFud's decision threshold is set to zero so that PolyFudBpp can analyze every putative PAS. (TODO: function to set the threshold)

This repository was created from [https://github.com/morikie/polar](https://github.com/morikie/polar). Since polar is too convoluted I separated PolyFud-BPP from it to clean it up and structure it a bit better.
