'''Calculates density of vessel types IT and PT.\n
   First calculates raw number of each vessel type\n
   by location IT or PT.'''

# Written by: Austin Schultz (aschultz37)
# Updated:    02/21/2022

import os
import numpy as np
import pandas as pd

class FileExtError(Exception):
    '''Catches error for wrong file type.'''
    pass 

def extract_file_tup(file_path):
    '''Returns a tuple of the file name and extension.'''
    file_basename = os.path.basename(file_path)
    return os.path.splitext(file_basename)

def dir_input():
    '''Returns the path to a directory.'''
    print('Enter the directory path:')
    dir_path = input()
    if(dir_path[-1] != '/'):
        dir_path = dir_path + '/'
    return dir_path

def read_file(file_path):
    '''Reads a file into a pandas dataframe and returns it.'''
    file_ext = extract_file_tup(file_path)[1]
    try:
        if file_ext == '.csv':
            return pd.read_csv(file_path)
        elif file_ext == '.xlsx':
            return pd.read_excel(file_path, index_col=0)
        else:
            raise FileExtError
    except FileNotFoundError:
        raise FileNotFoundError

def areafile_input():
    '''Gets the path to the file with IT/PT area for each ROI.'''
    print('Enter path to file with IT/PT areas:')
    areafile_path = input()
    return areafile_path

def output_dataframe(in_dir, file_path, df):
    '''Writes a pandas dataframe to a .csv file in 'output/' dir. New file\n
       name is original/input file name plus suffix _trimmed. Does not \n
       overwrite the original/input file.'''
    file_tup = extract_file_tup(file_path)
    file_ext = file_tup[1]
    outfile_path = 'output/' + in_dir + file_tup[0] + file_ext
    if file_ext == '.csv':
        df.to_csv(outfile_path, header=True, index=False)
    else:
        raise FileExtError    

def make_output_dir(in_dir):
    '''Checks if output directory exists. If not, makes the directory.\n
       The base directory name is 'output/in_dir/'.'''
    outfile_path = os.path.join('output/', in_dir)
    if os.path.exists(outfile_path) == False:
        os.makedirs(outfile_path)

def extract_sample_id(filename):
    '''Extracts the sample ID from file name.'''
    # Remove the "Feat_ROI_" and "_Final_trimmed"
    sample_num = filename[9:-14]
    return sample_num


def num_by_location(df, in_file):
    '''For an ROI, determines the number of each vessel type by location.\n
       Returns a dataframe listing counts of every vessel type for IT/PT.'''
    # Create new df listing: ROI name, # Type 1 (IT), ..., # Type 1 (PT), ...
    in_filename = extract_file_tup(in_file)[0]
    sample_num = extract_sample_id(in_filename)
    data_template = {'ROI': sample_num, '# Type 1 IT': 0, '# Type 2 IT': 0, 
                     '# Type 3 IT': 0, '# Type 4 IT': 0, '# Type 5 IT': 0,
                     '# Type 1 PT': 0, '# Type 2 PT': 0, '# Type 3 PT': 0, 
                     '# Type 4 PT': 0, '# Type 5 PT': 0}
    df_num = pd.DataFrame(data=data_template, index=[0])
    # Loop through df and increment appropriate col for each vessel
    df_t = df.T
    for col in df_t:
        vessel_type = str(df_t[col]['Vessel ID'])
        vessel_loc = df_t[col]['Vessel_Location']
        inc_string = '# Type ' + vessel_type + ' ' + vessel_loc[0].upper() + 'T'
        df_num.loc[:,(inc_string)] += 1
    return df_num

def find_it_area(sample_name, in_area):
    '''Finds the IT area for a sample and returns it.'''
    it_area = in_area.at[sample_name, 'IT Area']
    return it_area

def find_pt_area(sample_name, in_area):
    '''Finds the PT area for a sample and returns it.'''
    pt_area = in_area.at[sample_name, 'PT Area']
    return pt_area

def populate_areas(df, in_area):
    '''Creates columns that contain IT and PT area for each sample.'''
    df['IT Area'] = df['ROI'].map(lambda x: find_it_area(x, in_area))
    df['PT Area'] = df['ROI'].map(lambda x: find_pt_area(x, in_area))
    return df

def density_by_location(df):
    '''Using IT/PT areas from slides, calculate density of each vessel type\n
       for IT/PT. Returns a dataframe listing density of each vessel type\n
       for IT/PT.'''
    # Formula: (no. of vessels IT) / (total IT area) and same for PT
    # Append to df: Density Type 1 (IT), ..., Density Type 1 (PT), ...
    df['Rho Type 1 IT'] = df['# Type 1 IT'] / df['IT Area']
    df['Rho Type 2 IT'] = df['# Type 2 IT'] / df['IT Area']
    df['Rho Type 3 IT'] = df['# Type 3 IT'] / df['IT Area']
    df['Rho Type 4 IT'] = df['# Type 4 IT'] / df['IT Area']
    df['Rho Type 5 IT'] = df['# Type 5 IT'] / df['IT Area']
    df['Rho Type 1 PT'] = df['# Type 1 PT'] / df['PT Area']
    df['Rho Type 2 PT'] = df['# Type 2 PT'] / df['PT Area']
    df['Rho Type 3 PT'] = df['# Type 3 PT'] / df['PT Area']
    df['Rho Type 4 PT'] = df['# Type 4 PT'] / df['PT Area']
    df['Rho Type 5 PT'] = df['# Type 5 PT'] / df['PT Area']
    return df

def merge_roi(all_roi):
    '''Merges the calculations for all the ROI into one dataframe for\n
       easier output.'''
    if len(all_roi) == 0:
        print("Error: No ROI found.")
        return
    df_merged = pd.concat(all_roi, ignore_index=True)
    return df_merged


# Main program
in_dir = dir_input()
in_area_path = areafile_input()
in_area = read_file(in_area_path)
all_roi = list()
for in_file in os.listdir(in_dir):
    in_file_path = in_dir + in_file
    df = read_file(in_file_path)
    df_num = num_by_location(df, in_file)
    all_roi.append(df_num)
# merge all ROI into one dataframe
df_merged = merge_roi(all_roi)
# calculate density for each ROI
df_areas = populate_areas(df_merged, in_area)
df_density = density_by_location(df_areas)
# output the dataframe    
make_output_dir(in_dir)
output_dataframe(in_dir, 'vessel_density.csv', df_density)