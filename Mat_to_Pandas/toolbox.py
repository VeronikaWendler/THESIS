# Veronika Wendler
# 11.12.24

# Functions for data cleaning and wrangling and so on 
# 

##################
# import libraries
##################

import pandas as pd
import numpy as np
import re

#########################################################################################################################################################
# Functions for handling data conversion from .mat to pandas frames and for creating new columns and cleaning
#########################################################################################################################################################
# clean_column               -> unwraps and cleans categorical columns, nested lists
# parse_nested_list          -> parses nested string from info csv into obj, ensures formatting, this is used for the probability/image mapping
# rearrange_dataframe        -> rearranges the data so that op1 always contains the E option and op2 always contains the S option 
                                # even though they may have been positioned differently,
                                #other columns get updated according to this change (e.g. clicks etc)
#get_condition_garcia        -> creates condition mapping for the garcia replication


def clean_column(column, dtype='auto'):
    cleaned = []
    for x in column:
        while isinstance(x, (list, np.ndarray)) and len(x) > 0:
            x = x[0]
        if dtype == 'auto':
            dtype = 'numeric' if isinstance(x, (int, float, np.number)) else 'categorical'
        if dtype == 'numeric':
            try:
                cleaned.append(float(x))  
            except (ValueError, TypeError):
                cleaned.append(np.nan)  
        else:  
            cleaned.append(str(x))  
    return cleaned

def parse_nested_list(string):
    """ parses nested string from the info CSV into an object.
    Cleans up malformed structures and ensures proper formatting.
    """
    string = string.replace("array(", "").replace("dtype='<U5'", "").replace("dtype='<U6'", "").replace("dtype='<U4'", "").replace("dtype='<U12'", "").replace(")", "")
    string = re.sub(r",\s*,", ",", string)  
    string = re.sub(r"\[,\s*", "[", string) 
    string = re.sub(r",\s*\]", "]", string)  
    open_brackets = string.count("[")
    close_brackets = string.count("]")
    if open_brackets > close_brackets:
        string += "]" * (open_brackets - close_brackets)  
    elif close_brackets > open_brackets:
        string = "[" * (close_brackets - open_brackets) + string 
    return eval(string)  # eval

def clean_categorical_column(column):
        return [str(x[0]) if isinstance(x, (list, np.ndarray)) else str(x) for x in column]

def clean_numeric_column(column):
        """Cleans numeric columns by flattening nested arrays and converting to floats."""
        return [round(float(x[0]), 1) if isinstance(x, (list, np.ndarray)) else round(float(x), 1) for x in column]
    
# this is used in the EXP2...re files, it rearranges the position of particular values in particular columns
# op1 column only contains E and op2 only contains S options, other columns need to change their values accordingly
def rearrange_dataframe(df):
    swap_columns = ['cho', 'ev1', 'ev2', 'p1', 'p2', 'op1', 'op2', 'chose_right']
    for index, row in df.iterrows():
        if row['op1'] == 'S' and row['op2'] == 'E':
            # swap columns 'E' is in op1 and 'S' is in op2
            df.at[index, 'op1'], df.at[index, 'op2'] = row['op2'], row['op1']
            df.at[index, 'ev1'], df.at[index, 'ev2'] = row['ev2'], row['ev1']
            df.at[index, 'p1'], df.at[index, 'p2'] = row['p2'], row['p1']
            df.at[index, 'cho'] = 3 - row['cho']  
            df.at[index, 'chose_right'] = 1 - row['chose_right']  
    return df


#---------------------------------- Function specific to the garcia replication ----------------------------------
# creates the condition column
def get_condition_garcia(p1, p2):
    pair = sorted([p1, p2])
    if pair == [0.4, 0.6]:
        return 3  # 60/40
    elif pair == [0.3, 0.7]:
        return 2  # 70/30
    elif pair == [0.2, 0.8]:
        return 1  # 80/20
    elif pair == [0.1, 0.9]:
        return 0  # 90/10
    return -1
#------------------------------------ Function specific to the OV replication -------------------------------------------------------------------
# function that creates the condtion column
def get_OV_condition(p1, p2):
    pair = sorted([p1, p2])
    if pair == [0.1, 0.3]:
        return 3  # difficult
    elif pair == [0.2, 0.6]:
        return 2  # easy
    elif pair == [0.4, 0.8]:
        return 1  # easy
    elif pair == [0.7, 0.9]:
        return 0  # difficult
    return -1

#------------------------------------ Function specific to the Pilot replication -------------------------------------------------------------------
# cond col for Pilot Experiment
def get_pilot_condition(p1, p2):
    pair = sorted([p1, p2])
    if pair == [0.3750, 0.6250]:
        return 3  # 
    elif pair == [0.2500, 0.7500]:
        return 2  # 
    elif pair == [0.1250, 0.8750]:
        return 1  # 
    
    return -1