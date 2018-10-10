
import csv
import pandas as pd

# reading database
b_DB = pd.read_csv('biggDb.csv')

# reading exchange metabolites in liver
l_dt=pd.read_csv('livM_exMet.csv')

# reading metabolite annotation file
m_li = pd.read_csv('removed_compartment.csv')
m_li = pd.read_csv('metabolite_annotation_liver.csv')
m_a = pd.read_csv('rat-cobra.csv')

chem_ids = pd.merge(l_dt,m_li,how='left',on='cpd_id')

chem_ids.to_csv('chemIds.csv')

