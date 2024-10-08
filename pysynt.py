# PySynt Functions
# Louis Hague

import pandas as pd

def import_nucmer(file_path):
    coords_info = {
        "RStart": [],
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
                    coords_info["RStart"].append(int(line[0]))  
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

    # Calculate chromosome offsets
    fai_df['chrom_offset'] = fai_df['seq_len'].cumsum().shift(fill_value=0)

    return fai_df

import matplotlib.pyplot as plt
import matplotlib.patches as patches

def draw_segments_with_fill(RStart, Rend, Qstart, Qend, R_y=8, Q_y=4):
    # Create a figure and axis
    fig, ax = plt.subplots()
    
    # Define the corners of the filled area
    top_left = (RStart, R_y)
    top_right = (Rend, R_y)
    bottom_left = (Qstart, Q_y)
    bottom_right = (Qend, Q_y)
    
    # Create a polygon to fill the area between the lines
    filled_area = patches.Polygon([top_left, top_right, bottom_right, bottom_left], 
                                  closed=True, facecolor='lightblue', edgecolor='none', alpha=0.5)
    
    # Add the threads to the plot
    ax.add_patch(filled_area)
    
    ax.set_xlim(min(RStart, Qstart) - 1, max(Rend, Qend) + 1)
    ax.set_ylim(Q_y - 2, R_y + 2)
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_title('Filled Area Between Segments Without Top and Bottom Lines')
    
    # Add grid and reference lines
    plt.grid()
    plt.axhline(0, color='black', linewidth=0.5, ls='--')
    plt.axvline(0, color='black', linewidth=0.5, ls='--')
    
    # Display the plot
    plt.show()

# Example usage
draw_segments_with_fill(RStart=1, Rend=2, Qstart=2, Qend=3)