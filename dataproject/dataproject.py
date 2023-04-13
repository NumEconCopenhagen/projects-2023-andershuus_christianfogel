#These packages must be installed in order to create the maps. 
%pip install geopandas
%pip install mapclassify

import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets
from dstapi import DstApi

#We are getting data from two different tables from Statistics Denmark
church = DstApi('KM6')
inc = DstApi('INDKP132')

#We specify what we select from the dataset. We want the average family, and we want to see the total and 
#not specific income intervals.
params = {'table': 'indkf132',
 'format': 'BULK',
 'lang': 'en',
 'variables': [{'code': 'OMRÅDE', 'values': ['*']},
  {'code': 'ENHED', 'values': ['117']}, #Average income for families in the group (DKK)
  {'code': 'FAMTYP', 'values': ['*']},
  {'code': 'INDKINTB', 'values': ['99']}, #Total
  {'code': 'Tid', 'values': ['*']}]}

#We apply the dictionary created above to get our dataset.
inc_table = inc.get_data(params=params)

#We sort the data.
inc_table.sort_values(by=['OMRÅDE', 'TID', 'FAMTYP'], inplace=True)

#Removing non-used columns to simplify the data set
inc_table_000=inc_table.loc[:, ['OMRÅDE', 'TID', 'FAMTYP','INDHOLD']]
inc_table_000.head()










def keep_regs(df, regs):
    """ Example function. Keep only the subset regs of regions in data.

    Args:
        df (pd.DataFrame): pandas dataframe 

    Returns:
        df (pd.DataFrame): pandas dataframe

    """ 
    
    for r in regs:
        I = df.reg.str.contains(r)
        df = df.loc[I == False] # keep everything else
    
    return df