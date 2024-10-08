# PySynt Testing, Louis Hague

# Imports
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pysynt import import_genome, import_nucmer

# Change to data directory
os.chdir("C://Users//louis//syn_plot")

# Data
alignments = import_nucmer("sugar_cane//out.coords")
ref_data = import_genome("sugar_cane//Sspon.fna.fai")
query_data = import_genome("sugar_cane//CC_1940.fna.fai")
reference_scafs = ["Chr7C"]
query_scafs = ["CM034588.1","CM034589.1"]

# Subsetting Alignments for scaffolds of interest
alignments = alignments[(alignments['query'].isin(query_scafs)) & 
                                 (alignments['reference'].isin(reference_scafs))]

# Subset ref_data and query_data based on scafs 
ref_data = ref_data[(ref_data['seq_names'].isin(reference_scafs))]
query_data = query_data[(query_data['seq_names'].isin(query_scafs))]

# Initialise the plot
fig, ax = plt.subplots(figsize=(15, 8))

# First ref plot position
x_position = 0

# Plot the ref chromosomes as rectangles next to one another
for i, row in ref_data.iterrows():
    rect = patches.Rectangle((x_position, 8), row['seq_len'], 1, linewidth=1, edgecolor='black', facecolor='royalblue')
    ax.add_patch(rect)
    # Move the x pos for the start pos of the next rectangle
    x_position += row['seq_len']
    # Chromosome label
    ax.text(x_position - row['seq_len'] / 2, 8.5, row['seq_names'], ha='center')

# First query plot position
y_position = 0

# Plot the query chromosomes as rectangles next to one another
for i, row in query_data.iterrows():
    rect = patches.Rectangle((y_position, 4), row['seq_len'], 1, linewidth=1, edgecolor='black', facecolor='royalblue')
    ax.add_patch(rect)
    # Move the y pos for the start pos of the next rectangle
    y_position += row['seq_len']
    # Chromosome label
    ax.text(y_position - row['seq_len'] / 2, 4.5, row['seq_names'], ha='center')

# Plot details
ax.set_xlim(0, x_position)
ax.set_ylim(1, 10)
ax.set_title("Chromosome Blocks")
plt.show()