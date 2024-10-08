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