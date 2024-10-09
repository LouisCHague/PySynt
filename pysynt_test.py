# PySynt Testing, Louis Hague

# Imports
import os
from pysynt import import_genome, import_nucmer, plot_alignment_multi

# Change to data directory
os.chdir("C://Users//louis//syn_plot")

# Data
alignments = import_nucmer("melon//out.coords")
ref_data = import_genome("melon//Melon_351_.fasta.fai")
query_data = import_genome("melon//Melon_361_.fasta.fai")
reference_scafs = ["chr10"]
query_scafs = ["chr10"]
# Example usage
plot_alignment_multi(ref_data, query_data, reference_scafs, query_scafs, alignments, min_alignment_size=9000)