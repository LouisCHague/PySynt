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

# def reorder_chromosomes(alignments, ref_data, query_data):
#     # Step 1: Calculate average alignment length and create unique identifiers for reference-query pairs
#     alignments['average_length'] = alignments['Rlen'] + alignments['Qlen']
#     alignments['ref_query'] = alignments['reference'] + "_" + alignments['query']
#     alignment_median = alignments.groupby('ref_query')

#     print(ref_data)
#     print(query_data)
#     print(alignment_median)
    
#     return 

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


def create_chromosome(ax, genome_data, y_pos):
    '''Generates chromosome blocks'''

    # Initial position for the first chromosome
    x_position = 0

    # Iterate over each row in the genome_data DataFrame
    for i, row in genome_data.iterrows():
        # Create a rectangle for each chromosome
        rect = patches.Rectangle((x_position, y_pos), row['seq_len'], 1, linewidth=1, edgecolor='black', facecolor='grey')
        ax.add_patch(rect)
        
        # Position for the next chromosome
        x_position += row['seq_len']
        
        # Add chromosome label text
        ax.text(x_position - row['seq_len'] / 2, y_pos + 0.5, row['seq_names'], ha='center', va='center')
    return x_position





def plot_alignment_duo(ref_data, query_data, reference_scafs, query_scafs, alignments, min_alignment_size=5000):
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

    # Reorder the midpoints based on alignment
    #reorder_chromosomes(alignments,ref_data,query_data)
    # Offset changes where the alignments are positioned
    query_data['chrom_offset'] = query_data['seq_len'].cumsum().shift(fill_value=0)
    ref_data['chrom_offset'] = ref_data['seq_len'].cumsum().shift(fill_value=0)

    # Generate the chromosome blocks
    x_position = max(create_chromosome(ax, ref_data, 8), create_chromosome(ax, query_data, 4))


    # Plot alignment threads
    for i, row in alignments.iterrows():

        # Obtain offset for current chromosome (Ref + Query)
        off_ref = ref_data[ref_data['seq_names'] == row['reference']]
        off_ref = off_ref['chrom_offset']

        off_query = query_data[query_data['seq_names'] == row['query']]
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


######################################################
# Multiple Synteny Plotter
######################################################


def import_multi_coords(alignments,chromosomes):
    '''Imports coords files contained within a list and subsets them'''

    # Pointer (For alignments and chromosomes)
    a = 0
    b = 1

    # Dictionary to store dataframes on coords files
    coords_dict = {}

    # Subset the files
    for i in alignments:

        # Subset the coords files
        coords_file = import_nucmer(alignments[a])
        coords_file = coords_file[(coords_file['reference'].isin(chromosomes[a])) & 
                                 (coords_file['query'].isin(chromosomes[b]))]
        coords_dict[i] = coords_file 

        # Alter pointers for next chromosome
        a += 1
        b += 1
    return(coords_dict)

def import_multi_fai(genomes,chromosomes):
    '''Imports fai files contained within a list and subsets them'''

    # Dictionary to store dataframes on fai files
    fai_dict = {}

    # Subset the files
    for i in range(len(genomes)):
        genome_file = import_genome(genomes[i])
        genome_file = genome_file[(genome_file['seq_names'].isin(chromosomes[i]))]
        fai_dict[genomes[i]] = genome_file 
    return fai_dict




def plot_alignment_multi(genomes, alignments, chromosomes,min_alignment_size=5000):
    '''Plots alignment between multiple genomes'''

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Import bulk files
    alignment_df = import_multi_coords(alignments, chromosomes)
    fai_df = import_multi_fai(genomes,chromosomes)

    # Minimum alignment size
    for key, value in alignment_df.items():
        alignment_df[key] = value[
        (abs(value['Rend'] - value['Rstart']) >= min_alignment_size) & 
        (abs(value['Qend'] - value['Qstart']) >= min_alignment_size)
        ]
    
    # Index to pass through fai files
    n = 1
    # Order alignments
    for key, value in alignment_df.items():

        # Minimum alignment size
        alignment_df[key] = value[
        (abs(value['Rend'] - value['Rstart']) >= min_alignment_size) & 
        (abs(value['Qend'] - value['Qstart']) >= min_alignment_size)
        ]
        
        # fai_df information
        fai_keys = list(fai_df.keys())
        fai_values = list(fai_df.values())

        fai_df[fai_keys[n]]['chrom_offset'] = fai_values[n]['seq_len'].cumsum().shift(fill_value=0)

        # Iterate counter
        n += 1
    
    #for key,value in fai_df.items():
    #    print(value.head())

    # Generate chromosome blocks
    n = len(genomes) * 4
    # Need to make the x-axis as wide as the longest genome
    x_position = 0
    for key, value in fai_df.items():
        x_position = max(create_chromosome(ax, value, n), x_position)
        n -= 4
    
    # Need to offset the chromosomes for the first genome
    fai_df[fai_keys[0]]['chrom_offset'] = fai_values[0]['seq_len'].cumsum().shift(fill_value=0)
    fai_keys = list(fai_df.keys())
    fai_values = list(fai_df.values())

    # print(fai_values)

    # Plot alignment threads
    n = 0
    y_axis_coord = len(genomes) * 4

    for key, value in alignment_df.items():
        for i, row in value.iterrows():

            #print(fai_values)
            #print(value)

            # Obtain offset for current chromosome (Ref + Query)
            if n == 0:
                off_ref = fai_values[n][fai_values[n]['seq_names'] == row['reference']]
                #print(f'off_ref: {off_ref}')
                off_ref = off_ref['chrom_offset']

                # Only those which have been altered will have query_order value
                off_query = fai_values[n + 1][fai_values[n + 1]['seq_names'] == row['query']]
                #print(f'off_query: {off_query}')
                off_query = off_query['chrom_offset']         
            else:
                # The original genome is fixed, the rest are reordered
                #print(row)
                #print(fai_values[n]['query_order'])
                off_ref = fai_values[n][fai_values[n]['seq_names'] == row['reference']]
                #print(f'off_ref: {off_ref}')
                off_ref = off_ref['chrom_offset']   

                #print(row)
                off_query = fai_values[n + 1][fai_values[n + 1]['seq_names'] == row['query']]
                #print(f'off_query: {off_query}')
                off_query = off_query['chrom_offset']   
        
            # Add offset value to the thread alignment coordinate
            Rstart = row['Rstart'] + int(off_ref.iloc[0])
            Rend = row['Rend'] + int(off_ref.iloc[0])
            Qstart = row['Qstart'] + int(off_query.iloc[0])
            Qend = row['Qend'] + int(off_query.iloc[0])

            # Generate thread polygon
            if n == 0:
                thread = create_thread(Rstart, Rend, Qstart, Qend, R_y=y_axis_coord, Q_y=y_axis_coord - 3)
                ax.add_patch(thread)
            else:
                thread = create_thread(Rstart, Rend, Qstart, Qend, R_y=y_axis_coord, Q_y=y_axis_coord - 3)
                ax.add_patch(thread)

        # Change the genome pair we are iterating over
        # each time the alignment changes
        n += 1
        y_axis_coord -= 4
    
    # Graph Modifiers 
    ax.set_xlim(0, x_position + 100)
    ax.set_ylim(3, len(genomes) * 4 + 2)
    ax.set_title("Chromosome Blocks")
    ax.set_xlabel("Base Pair Length")
    #plt.yticks([])
    ax.ticklabel_format(useOffset=False, style='plain')
    plt.grid()
    plt.axhline(0, color='black', linewidth=0.5, ls='--')
    plt.axvline(0, color='black', linewidth=0.5, ls='--')
    # Display the plot
    plt.show()