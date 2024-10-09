# PySynt Testing, Louis Hague

# Imports
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pysynt import import_genome, import_nucmer, calculate_midpoints, reorder_query_chromosomes, create_thread

# Change to data directory
os.chdir("C://Users//louis//syn_plot")

# Data
alignments = import_nucmer("sugar_cane//out.coords")
ref_data = import_genome("sugar_cane//Sspon.fna.fai")
query_data = import_genome("sugar_cane//CC_1940.fna.fai")
reference_scafs = ["Chr7C"]
query_scafs = ["CM034588.1","CM034589.1"]

###################################################################################
# Need to make a subsetting function for ease or append it to the plotting function 
###################################################################################

# Subsetting Alignments for scaffolds of interest
alignments = alignments[(alignments['query'].isin(query_scafs)) & 
                                 (alignments['reference'].isin(reference_scafs))]
# Subset ref_data and query_data based on scafs 
ref_data = ref_data[(ref_data['seq_names'].isin(reference_scafs))]
query_data = query_data[(query_data['seq_names'].isin(query_scafs))]

###################
# Plotting function
###################

def plot_alignment_multi(ref_data, query_data, alignments):
    '''Synteny plotting function'''

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Calculate midpoints, reorders query_data and alignments
    sorted_queries, alignments = calculate_midpoints(alignments)
    query_data = reorder_query_chromosomes(query_data, sorted_queries)

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

    
    # Plot alignment threads
    for i, row in alignments.iterrows():
        # Need to offset the threads by the previous chromosome 
        # MIGHT HAVE TO DO THIS FOR REFERENCE CHROMOSOMES TOO! 
        off_row = query_data[query_data['query_order'] == row['query']]
        offset = off_row['chrom_offset']
        # print(offset + 2)
        
        # Add offset to the threads
        Rstart = row['Rstart'] + int(offset.iloc[0])
        Rend = row['Rend'] + int(offset.iloc[0])
        Qstart = row['Qstart'] + int(offset.iloc[0])
        Qend = row['Qend'] + int(offset.iloc[0])

        # Generate thread polygon
        thread = create_thread(Rstart, Rend, Qstart, Qend)
        ax.add_patch(thread)
    
    ax.set_xlim(0, x_position)
    ax.set_ylim(1, 10)
    ax.set_title("Chromosome Blocks")
    plt.grid()
    plt.axhline(0, color='black', linewidth=0.5, ls='--')
    plt.axvline(0, color='black', linewidth=0.5, ls='--')
    # Display the plot
    plt.show()

# Example usage
plot_alignment_multi(ref_data, query_data, alignments.head(500))