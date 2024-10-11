# PySynt Testing, Louis Hague

# Imports
import os
from pysynt import import_genome, import_nucmer, plot_alignment_duo, plot_alignment_multi

# Change to data directory
os.chdir("C://Users//louis//syn_plot")

# Data
alignments = import_nucmer("melon//out.coords")
ref_data = import_genome("melon//Melon_351_.fasta.fai")
#print(ref_data)
query_data = import_genome("melon//Melon_361_.fasta.fai")
reference_scafs = ["chr1","chr10"]
query_scafs = ["chr01","chr10"]
# Example usage
#plot_alignment_duo(ref_data, query_data, reference_scafs, query_scafs, alignments, min_alignment_size=10000)


alignments = import_nucmer("cancer//out.coords")
ref_data = import_genome("cancer//GRCh38.fna.fai")
#print(ref_data)
query_data = import_genome("cancer//HCC.fna.fai")
reference_scafs = ["NC_000006.12"]
query_scafs = ["scaffold_5","scaffold_30","scaffold_41"]
# Example usage
plot_alignment_duo(ref_data, query_data, reference_scafs, query_scafs, alignments, min_alignment_size=10000)

# List of coords files
#alignments = ["Yeast//alignments//out_146_211.coords", "Yeast//alignments//out_211_235.coords"]
# Genomes files 
#genomes = ["Yeast//genomes//GCA_000146045.2_R64_genomic.fna.fai",
#             "Yeast//genomes//GCA_021172205.1_ASM2117220v1_genomic.fna.fai",
#             "Yeast//genomes//GCA_023508825.1_ASM2350882v1_genomic.fna.fai"]
# Chromosomes of interest (First genome scaffs)...
#chromosomes = [["BK006935.2"], ["CP089100.1"], ["CP096554.1"]]

# Example usage
#plot_alignment_multi(genomes, alignments, chromosomes)