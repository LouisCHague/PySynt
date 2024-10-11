# PySynt Testing
# Louis Hague

# Imports
import os
from pysynt import import_genome, import_nucmer, plot_alignment_duo, plot_alignment_multi

# Change to data directory
os.chdir("C://Users//louis//syn_plot")

# Data
alignments = import_nucmer("melon//out.coords")
ref_data = import_genome("melon//Melon_351_.fasta.fai")
query_data = import_genome("melon//Melon_361_.fasta.fai")
# Chromosomes of interest
reference_scafs = ["chr1","chr10"]
query_scafs = ["chr01","chr10"]
# Example usage
plot_alignment_duo(ref_data, query_data, reference_scafs, query_scafs, alignments, min_alignment_size=10000)

# List coords alignment files
alignments = ["Yeast_2//alignments//out_146_211.coords", "Yeast_2//alignments//out_211_235.coords"]
# Genome files 
genomes = ["Yeast_2//genomes//GCF_000146_sample.fna.fai",
             "Yeast_2//genomes//GCA_0211_sample.fna.fai",
             "Yeast_2//genomes//GCA_0235_sample.fna.fai"]
# Chromosomes of interest
chromosomes = [["NC_001133.9", "NC_001134.8"], ["CP089100.1", "CP089101.1"], ["CP096554.1", "CP096553.1"]]
# Example usage
plot_alignment_multi(genomes, alignments, chromosomes)