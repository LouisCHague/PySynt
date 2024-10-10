# PySynt Testing, Louis Hague

# Imports
import os
from pysynt import import_genome, import_nucmer, plot_alignment_duo, plot_alignment_multi, import_multi_fai

# Change to data directory
os.chdir("C://Users//louis//syn_plot")

# Data
alignments = import_nucmer("melon//out.coords")
ref_data = import_genome("melon//Melon_351_.fasta.fai")
#print(ref_data)
query_data = import_genome("melon//Melon_361_.fasta.fai")
reference_scafs = ["chr10"]
query_scafs = ["chr10"]
# Example usage
#plot_alignment_duo(ref_data, query_data, reference_scafs, query_scafs, alignments, min_alignment_size=9000)



# I could make it so as many genomes as the user puts in are plotted
# It would mean that we need to organise the (n+1) of whatever genome is
# put above it similar to how I organised the query in the current example
# I would also need to set a way to alter the y-axis of the two genomes 
# and there alignments based on how many are set. A block of chromosomes
# takes up 1 on the y-axis as well as the space between them. 
# Y-axis space:
# 2 blocks, 1 space = 3
# 3 blocks, 2 spaces = 5
# 4 blocks, 3 spaces = 7

# List of coords files
alignments = ["Yeast//alignments//out_146_211.coords", "Yeast//alignments//out_211_235.coords"]
# Genomes files 
genomes = ["Yeast//genomes//GCA_000146045.2_R64_genomic.fna.fai",
             "Yeast//genomes//GCA_021172205.1_ASM2117220v1_genomic.fna.fai",
             "Yeast//genomes//GCA_023508825.1_ASM2350882v1_genomic.fna.fai"]
# Chromosomes of interest (First genome scaffs)...
chromosomes = [["BK006935.2"], ["CP089100.1"], ["CP096554.1"]]

# Example usage
plot_alignment_multi(genomes, alignments, chromosomes)