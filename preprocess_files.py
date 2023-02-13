'''Filters vessels based on area.\n
   Then calculates vessel type ID based on\n
   marker thresholds. Input and output are .csv.'''

# Written by: Austin Schultz (aschultz37)
# Updated:    02/13/2022

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
    '''Reads a .csv file into a pandas dataframe and returns it.'''
    file_ext = extract_file_tup(file_path)[1]
    try:
        if file_ext == '.csv':
            return pd.read_csv(file_path)
        else:
            raise FileExtError
    except FileNotFoundError:
        raise FileNotFoundError

def output_dataframe(in_dir, file_path, df):
    '''Writes a pandas dataframe to a .csv file in 'output/' dir. New file\n
       name is original/input file name plus suffix _trimmed. Does not \n
       overwrite the original/input file.'''
    file_tup = extract_file_tup(file_path)
    file_ext = file_tup[1]
    outfile_path = 'output/' + in_dir + file_tup[0] + '_trimmed' + file_ext
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

def do_trim(df):
    '''Removes items from the dataframe based on given criteria.\n
       Returns a copy of the original dataframe (with items removed).'''
    area_threshold = 600
    df_t = df.T
    for col in df_t.columns:
        vessel_area = df_t[col][0]
        #print(vessel_area)
        if int(vessel_area) < area_threshold:
            #print("pop " + str(col))
            df_t.pop(col)
            col = col-1
    return df_t.T

def calc_vessel_type(df):
    '''Calculates vessel ID based on T/F for marker positivity. Marker \n
       positivity determined based on pre-set thresholds. Adds columns \n
       to dataframe and returns.'''
    d240_threshold = 0.1                # if >, Lymphatic; if <, Blood
    lyve1_threshold = 0.7               # Lymphatic: if >, 2; if <, 1
    asma_threshold = 0.45               # Blood: if >, add 1; if <, add 2
    cd34_threshold = 2                  # Blood: if >, add 3; if <, add 2
    aqp1_threshold = 0.5                # Blood: if >, 'True'; if <, 'False'
    df['Vessel Type'] = df['Int_D240'].map(lambda x: 
                                           'Lymphatic' if x > d240_threshold
                                           else 'Blood')
    df['LYVE1+/-'] = df['Int_LYVE1'].map(lambda x:
                                         2 if x > lyve1_threshold
                                         else 1)
    df.loc[df['Vessel Type'] == 'Blood', 'LYVE1+/-'] = 0
    df['aSMA+/-'] = df['Int_aSMA'].map(lambda x:
                                       1 if x > asma_threshold
                                       else 2)
    df.loc[df['Vessel Type'] == 'Lymphatic', 'aSMA+/-'] = 0
    df['CD34+/-'] = df['Int_CD34'].map(lambda x:
                                       3 if x > cd34_threshold
                                       else 2)
    df.loc[df['Vessel Type'] == 'Lymphatic', 'CD34+/-'] = 0
    df['Vessel ID'] = df.loc[:, 'LYVE1+/-':'CD34+/-'].sum(1)
    df['AQP1'] = df['Int_AQP1'].map(lambda x:
                                    'True' if x > aqp1_threshold
                                    else 'False')
    df.loc[df['Vessel Type'] == 'Lymphatic', 'AQP1'] = 0
    return df

in_dir = dir_input()
for in_file in os.listdir(in_dir):
    in_file_path = in_dir + in_file
    df = read_file(in_file_path)
    df_trimmed = do_trim(df)
    df_trimmed_calc = calc_vessel_type(df_trimmed)
    make_output_dir(in_dir)
    output_dataframe(in_dir, in_file, df_trimmed_calc)
