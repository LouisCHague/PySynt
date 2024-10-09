# PySynt Functions
# Louis Hague

# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def import_nucmer(file_path):
    '''Imports a pd.df of a nucmer coords file'''

    # Dict. to hold coords info
    coords_info = {
        "Rstart": [],
        "Rend": [],
        "Qstart": [],
        "Qend": [],
        "Rlen": [],
        "Qlen": [],
        "identity": [],
        "reference": [],
        "query": []
    }
    
    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip().split("\t")
                
                if len(line) == 9:
                    coords_info["Rstart"].append(int(line[0]))  
                    coords_info["Rend"].append(int(line[1]))    
                    coords_info["Qstart"].append(int(line[2]))   
                    coords_info["Qend"].append(int(line[3]))    
                    coords_info["Rlen"].append(int(line[4])) 
                    coords_info["Qlen"].append(int(line[5]))   
                    coords_info["identity"].append(float(line[6])) 
                    coords_info["reference"].append(line[7])   
                    coords_info["query"].append(line[8])    

        coords_df = pd.DataFrame(coords_info)
        
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        return pd.DataFrame() 
    except Exception as e:
        print(f"An error occurred: {e}")
        return pd.DataFrame()

    return coords_df

def get_contig_lengths(fai_file):
    '''Gets the lengths of the contigs'''

    seq_len = {}
    with open(fai_file, 'r') as file:
        for line in file:
            fields = line.strip().split("\t")
            seq_len[fields[0]] = int(fields[1])
    return seq_len

def import_genome(fai_file):
    '''Imports information from fai files'''

    seq_len = get_contig_lengths(fai_file)

    fai_df = pd.DataFrame({
        'seq_names': seq_len.keys(),
        'seq_len': seq_len.values(),
        # Assume all orientations are '+'
        'seq_ori': ['+'] * len(seq_len),  
    })

    return fai_df

######################################################################
# Reorders the alignments and query_data so the plot is easier to read
######################################################################

def calculate_midpoints(alignments):
    '''Calculate midpoint for each alignment, generates ordered list of queries
    AND reorders the alignments dataframe based on ordered query chromosomes'''

    # Calculate the midpoint for each alignment
    alignments['ref_midpoint'] = (alignments['Rstart'] + alignments['Rend']) / 2
    
    # Group by query chromosome, then calculate median midpoint for each group
    query_midpoints = alignments.groupby('query')['ref_midpoint'].median()
    
    # Sort query chromosomes based on the calculated midpoints
    sorted_queries = query_midpoints.sort_values().index.tolist()

    # I want to convert the alignments so they follow the new order of the chromosomes
    alignments['query'] = pd.Categorical(alignments['query'], categories=sorted_queries, ordered=True)
    alignments = alignments.sort_values(by='query').reset_index(drop=True)
    
    return sorted_queries, alignments

def reorder_query_chromosomes(query_data, sorted_queries):
    '''Reorders query_data based on new order'''

    query_data['query_order'] = pd.Categorical(query_data['seq_names'], categories=sorted_queries, ordered=True)
    query_data = query_data.sort_values('query_order')

    # Need to create an offset to offset the alignment strands in the graph
    query_data['chrom_offset'] = query_data['seq_len'].cumsum().shift(fill_value=0)
    
    return query_data

#############################
# Plots the alignment threads
#############################

def create_thread(Rstart, Rend, Qstart, Qend, R_y=8, Q_y=5):
    '''Generates the synteny thread'''
    
    # Define the corners of the filled area
    top_left = (Rstart, R_y)
    top_right = (Rend, R_y)
    bottom_left = (Qstart, Q_y)
    bottom_right = (Qend, Q_y)

    # Different colours for reverse sequences
    if Rstart > Rend or Qstart > Qend:
        filled_area = patches.Polygon([top_left, top_right, bottom_right, bottom_left], 
                                  closed=True, facecolor='red', edgecolor='none', alpha=0.5)
    else:
        filled_area = patches.Polygon([top_left, top_right, bottom_right, bottom_left], 
                                  closed=True, facecolor='lightblue', edgecolor='none', alpha=0.5)
    
    return filled_area

#fig, ax = plt.subplots()
#polygon = create_thread(Rstart=1,Rend=3,Qstart=1,Qend=4)
#ax.add_patch(polygon)
#ax.set_title("Chromosome Blocks")
#plt.grid()
#ax.set_xlim(0, 10)
#ax.set_ylim(1, 10)
#plt.axhline(0, color='black', linewidth=0.5, ls='--')
#plt.axvline(0, color='black', linewidth=0.5, ls='--')
# Display the plot
#plt.show()

def plot_alignment_multi(ref_data, query_data, reference_scafs, query_scafs, alignments, min_alignment_size=5000):
    '''Synteny plotting function'''

    # Subsetting Alignments for scaffolds of interest
    alignments = alignments[(alignments['query'].isin(query_scafs)) & 
                                 (alignments['reference'].isin(reference_scafs))]
    # Subset ref_data and query_data based on scafs 
    ref_data = ref_data[(ref_data['seq_names'].isin(reference_scafs))]
    query_data = query_data[(query_data['seq_names'].isin(query_scafs))]

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
        rect = patches.Rectangle((x_position, 8), row['seq_len'], 1, linewidth=1, edgecolor='black', facecolor='grey')
        ax.add_patch(rect)
        # Move the x pos for the start pos of the next rectangle
        x_position += row['seq_len']
        # Chromosome label
        ax.text(x_position - row['seq_len'] / 2, 8.5, row['seq_names'], ha='center')

    # First query plot position
    y_position = 0

    # Plot the query chromosomes as rectangles next to one another
    for i, row in query_data.iterrows():
        rect = patches.Rectangle((y_position, 4), row['seq_len'], 1, linewidth=1, edgecolor='black', facecolor='grey')
        ax.add_patch(rect)
        # Move the y pos for the start pos of the next rectangle
        y_position += row['seq_len']
        # Chromosome label
        ax.text(y_position - row['seq_len'] / 2, 4.5, row['seq_names'], ha='center')

    # Need to calculate offset for reference chromosomes
    # This step was done in reorder_query_chromosomes for query chromosomes
    ref_data['chrom_offset'] = ref_data['seq_len'].cumsum().shift(fill_value=0)

    # Plot alignment threads
    for i, row in alignments.iterrows():

        # Obtain offset for current chromosome (Ref + Query)
        off_ref = ref_data[ref_data['seq_names'] == row['reference']]
        off_ref = off_ref['chrom_offset']

        off_query = query_data[query_data['query_order'] == row['query']]
        off_query = off_query['chrom_offset']
        # print(offset + 2)
        
        # Add offset value to the thread alignment coordinate
        Rstart = row['Rstart'] + int(off_ref.iloc[0])
        Rend = row['Rend'] + int(off_ref.iloc[0])
        Qstart = row['Qstart'] + int(off_query.iloc[0])
        Qend = row['Qend'] + int(off_query.iloc[0])

        # Generate thread polygon
        thread = create_thread(Rstart, Rend, Qstart, Qend)
        ax.add_patch(thread)
    
    # Graph Modifiers 
    ax.set_xlim(0, x_position + 100)
    ax.set_ylim(3, 10)
    ax.set_title("Chromosome Blocks")
    ax.set_xlabel("Base Pair Length")
    plt.yticks([])
    ax.ticklabel_format(useOffset=False, style='plain')
    plt.grid()
    plt.axhline(0, color='black', linewidth=0.5, ls='--')
    plt.axvline(0, color='black', linewidth=0.5, ls='--')
    # Display the plot
    plt.show()