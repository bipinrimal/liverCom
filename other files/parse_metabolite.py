#Reading Json file
import json
import re
import csv
import pandas as pd

def find_word(text,search):
    result = re.findall('\\b'+search+'\\',text, flags=re.IGNORECASE)
    if len(result)>0:
        return True
    else:
        return False

with open ('universal_model.json','r') as um:
    obj = json.load(um)

df=pd.DataFrame()
for metabolite in obj['metabolites']:
    annotation = metabolite['annotation']
    for text in annotation:
        if "CHEBI" in text:
            chebi = text[1].replace('http://identifiers.org/chebi/CHEBI:','')
        if "KEGG Compound" in text:
            kegg = text[1].replace('http://identifiers.org/kegg.compound/','')
    values=[[metabolite['id'],chebi,kegg,metabolite['name']]]
    df = df.append(values)


df.to_csv('biggDb.csv', sep=',')
    #print(metabolite['id'],metabolite['name'],chebi,kegg)


