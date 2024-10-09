# PySynt: Python functions for exploring synteny using whole genome alignments

## Data Pre-processing

Index your reference and query Fasta files with [SAMtools faidx](https://www.htslib.org/doc/samtools-faidx.html):
```samtools faidx ref.fasta```

Align your reference and query Fasta files with [MUMmer](https://mummer4.github.io/):
```nucmer [options] <reference file> <query file>```
```show-coords -T -H -r ref_qry.delta > ref_qry.coords```

## Inspiration
This project was inspired by the far more sophisticated [asynt functions](https://github.com/simonhmartin/asynt/tree/master) project.
