'''Calculates density of vessel types IT and PT.\n
   First calculates raw number of each vessel type\n
   by location IT or PT.\n
   Expects file names in format ##-###_[<ROI>].csv!'''

# Written by: Austin Schultz (aschultz37)
# Updated:    03/29/2023

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
        df.to_csv(outfile_path, header=True, index=True)
    else:
        raise FileExtError    

def make_output_dir(in_dir):
    '''Checks if output directory exists. If not, makes the directory.\n
       The base directory name is 'output/in_dir/'.'''
    outfile_path = os.path.join('output/', in_dir)
    if os.path.exists(outfile_path) == False:
        os.makedirs(outfile_path)

def extract_sample_id(filename):
    '''Extracts the sample ID and ROI from file name. Returns as tuple.\n
       Expects file name in format ##-###_[<ROI>].csv'''
    sample_num = filename[0:6]
    roi = filename[7:]
    return (sample_num, roi)


def num_by_location(df, in_file):
    '''For an ROI, determines the number of each vessel type by location.\n
       Returns a dataframe listing counts of every vessel type for IT/PT.'''
    # Create new df listing: ROI name, # Type 1 (IT), ..., # Type 1 (PT), ...
    in_filename = extract_file_tup(in_file)[0]
    sample_num = extract_sample_id(in_filename)
    data_template = {'Pt Number': sample_num[0],
                     'Sample ID': sample_num[0] + '_' + sample_num[1],
                     '# Type 1 IT': 0, '# Type 2 IT': 0, 
                     '# Type 3 IT': 0, '# Type 4 IT': 0, '# Type 5 IT': 0,
                     '# Type 1 PT': 0, '# Type 2 PT': 0, 
                     '# Type 3 PT': 0, '# Type 4 PT': 0, '# Type 5 PT': 0}
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
    try:
        it_area = in_area.at[sample_name, 'IT Area']
    except Exception:
        it_area = -128 # abs(-128) is >> other areas, will be obviously wrong
    return it_area

def find_pt_area(sample_name, in_area):
    '''Finds the PT area for a sample and returns it.'''
    try:
        pt_area = in_area.at[sample_name, 'PT Area']
    except Exception:
        pt_area = -128 # abs(-128) is >> other areas, will be obviously wrong
    return pt_area

def populate_areas(df, in_area):
    '''Creates columns that contain IT and PT area for each sample.'''
    df['IT Area'] = df['Sample ID'].map(lambda x: find_it_area(x, in_area))
    df['PT Area'] = df['Sample ID'].map(lambda x: find_pt_area(x, in_area))
    return df

def density_by_location(df):
    '''Using IT/PT areas from slides, calculate density of each vessel type\n
       for IT/PT. Returns a dataframe listing density of each vessel type\n
       for IT/PT.'''
    # Formula: (no. of vessels IT) / (total IT area) and same for PT
    # Append to df: Density Type 1 (IT), ..., Density Type 1 (PT), ...
    num_vessel_types = 5
    for i in range(1, num_vessel_types+1, 1):
        idx = str(i)
        df['Rho Type ' + idx + ' IT'] = (df['# Type ' + idx + ' IT'] / 
                                         df['IT Area'])
    for i in range(1, num_vessel_types+1, 1):
        idx = str(i)
        df['Rho Type ' + idx + ' PT'] = (df['# Type ' + idx + ' PT'] / 
                                         df['PT Area'])
    return df

def merge_roi(all_roi):
    '''Merges the calculations for all the ROI into one dataframe for\n
       easier output.'''
    if len(all_roi) == 0:
        print("Error: No ROI found.")
        return
    df_merged = pd.concat(all_roi, ignore_index=True)
    return df_merged

def merge_patients(df):
    '''Merges all ROI into one entry per patient by adding columns.'''
    df_merged_pt = df.groupby(by='Pt Number').sum(numeric_only=False)
    return df_merged_pt

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
df_merged_roi = merge_roi(all_roi)
# find the areas for each ROI
df_areas = populate_areas(df_merged_roi, in_area)
# merge all entries by patient
df_merged_pt = merge_patients(df_areas)
# calculate vessel densities for each patient
df_density = density_by_location(df_merged_pt)
# output the dataframe    
make_output_dir(in_dir)
output_dataframe(in_dir, 'vessel_density.csv', df_density)