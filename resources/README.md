# Reference resources

miR3End does not redistribute organism reference sequences or annotations.
Prepare the following resources for the organism and annotation release used
in your analysis, and pass their paths on the command line (or through the
equivalent `MIR3END_*` environment variables):

- miRNA hairpin Bowtie 1 index prefix: `--mir-hairpin`
- tRNA/snoRNA Bowtie 1 index prefix: `--trsno`
- genome Bowtie 1 index prefix: `--genome-index`
- genome FASTA: `--genome-fasta`
- GFF3 with `gene`/`ncRNA_gene` biotype attributes: `--rnatype-annotation`
- miRNA start/length table: `--meta-file`
- tail-base mechanism and sequence-merge tables: `--mechanism-file` and
  `--sequence-merge-file`

Record the source URL, release number and checksum of every reference file in
the analysis metadata. The small synthetic files under `test/data/` are only
for the included smoke test and are not biological reference data.
