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
alignments = import_nucmer("giardia//out.coords")
ref_data = import_genome("giardia//WBC6_4_5.fna.fai")
query_data = import_genome("giardia//BE2_4_5.fna.fai")
reference_scafs = ["NC_051860.1", "NC_051859.1"]
query_scafs = ["CP110919.1", "CP110920.1"]

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

#print(ref_data.head())
#print(query_data.head())
#print(alignments.head())

def plot_alignment_multi(ref_data, query_data, alignments, min_alignment_size):
    '''Synteny plotting function'''

    # Minimum alignment size
    alignments = alignments[
    (abs(alignments['Rend'] - alignments['Rstart']) >= min_alignment_size) & 
    (abs(alignments['Qend'] - alignments['Qstart']) >= min_alignment_size)
    ]

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

    # Need to offset reference alignments
    # This step was done in reorder_query_chromosomes for query 
    ref_data['chrom_offset'] = ref_data['seq_len'].cumsum().shift(fill_value=0)

    # Plot alignment threads
    for i, row in alignments.iterrows():
        off_ref = ref_data[ref_data['seq_names'] == row['reference']]
        off_ref = off_ref['chrom_offset']

        off_query = query_data[query_data['query_order'] == row['query']]
        off_query = off_query['chrom_offset']
        # print(offset + 2)
        
        # Add offset to the threads
        Rstart = row['Rstart'] + int(off_ref.iloc[0])
        Rend = row['Rend'] + int(off_ref.iloc[0])
        Qstart = row['Qstart'] + int(off_query.iloc[0])
        Qend = row['Qend'] + int(off_query.iloc[0])

        # Generate thread polygon
        thread = create_thread(Rstart, Rend, Qstart, Qend)
        ax.add_patch(thread)
    
    ax.set_xlim(0, x_position + 100)
    ax.set_ylim(1, 10)
    ax.set_title("Chromosome Blocks")
    plt.grid()
    plt.axhline(0, color='black', linewidth=0.5, ls='--')
    plt.axvline(0, color='black', linewidth=0.5, ls='--')
    # Display the plot
    plt.show()

# Example usage
plot_alignment_multi(ref_data, query_data, alignments, 9000)