import csv
import pandas as pd
import re


l_dt = pd.read_csv('livM_exMet.csv', header=None)
l_dt= l_dt.drop(l_dt.columns[1],axis = 1)
l_dt.columns =['Abbreviation']

ma_liv = pd.read_csv('metabolite_annotations.csv')
ma_liv=ma_liv.drop(columns=['Neutral formula','Charged formula', 'Charge','InChI string', 'SMILES'],axis=1)



chemIDs = pd.merge(l_dt,ma_liv,how='left',on='Abbreviation')
print(chemIDs)

kegg_IDs=[]
count = 0
for index, row in chemIDs.iterrows():
    if row['ChEBI ID']!='nan':
        if row['KEGG ID']=='nan':
           kegg_IDs.append(str('%s%s'%(row['ChEBI ID'],row['Abbreviation'][-3:])))

        if row ['KEGG ID']!='nan':
            kegg_IDs.append(str('%s%s' % (row['KEGG ID'], row['Abbreviation'][-3:])))
        else:
            kegg_IDs.append(str())
    print(index)
print(kegg_IDs)




