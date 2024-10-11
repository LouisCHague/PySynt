# PySynt: Python functions for exploring synteny using whole genome alignments

## Data Pre-processing

Index your reference and query Fasta files with [SAMtools faidx](https://www.htslib.org/doc/samtools-faidx.html)
* ```samtools faidx ref.fasta```

Align your reference and query Fasta files with [MUMmer](https://mummer4.github.io/)
* ```nucmer [options] <reference file> <query file>```
* ```show-coords -T -H -r ref_qry.delta > ref_qry.coords```

## Visualisation ([Data](https://www.nature.com/articles/s41598-018-26416-2))

* plot_alignment_duo(ref_data, query_data, reference_scafs, query_scafs, alignments, min_alignment_size=10000)
![Figure_1](https://github.com/user-attachments/assets/2686c927-593d-45d7-a51b-4e9ef3085451)

* plot_alignment_multi(genomes, alignments, chromosomes)
![multiple_alignment_plot](https://github.com/user-attachments/assets/208bb872-2cc3-4fd2-b4b5-2e446c4eb28d)

## Inspiration
This project was inspired by the far more sophisticated [asynt functions](https://github.com/simonhmartin/asynt/tree/master) project.
