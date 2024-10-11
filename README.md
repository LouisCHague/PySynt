# PySynt: Python functions for visualising whole genome alignments

## Visualising the alignment of two genomes
```plot_alignment_duo(ref_data, query_data, reference_scafs, query_scafs, alignments, min_alignment_size=10000)```

### Data Pre-processing

Index your reference and query Fasta files with [SAMtools faidx](https://www.htslib.org/doc/samtools-faidx.html)
* ```samtools faidx ref.fasta```

Align your reference and query Fasta files with [MUMmer](https://mummer4.github.io/)
* ```nucmer [options] <reference file> <query file>```
* ```show-coords -T -H -r ref_qry.delta > ref_qry.coords```

![Figure_1](https://github.com/user-attachments/assets/2686c927-593d-45d7-a51b-4e9ef3085451)

## Visualising the alignment of multiple genomes
```plot_alignment_multi(genomes, alignments, chromosomes)```

### Data Pre-processing

Index your Fasta files with [SAMtools faidx](https://www.htslib.org/doc/samtools-faidx.html)
* ```samtools faidx ref.fasta```

Align your Fasta files (n, n+1) with [MUMmer](https://mummer4.github.io/)
* ```nucmer [options] <reference file> <query file>```
* ```show-coords -T -H -r ref_qry.delta > ref_qry.coords```

![multiple_alignment_plot](https://github.com/user-attachments/assets/208bb872-2cc3-4fd2-b4b5-2e446c4eb28d)

## Inspiration
This project was inspired by the far more sophisticated [asynt functions](https://github.com/simonhmartin/asynt/tree/master) project.
