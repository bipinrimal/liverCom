# liverCom
microbiome-liver model


**************
### Hi Joshua,

I tried with the kegg identifier and making all of the metabolites similar. Parsing of the biggDb was easy and I created a easier csv file so that extraction will be easier next time. \
I am only working on the exchange metabolites. Making them uniform across all the models should be easy.\
But the liver model was difficult to work with. The number of metabolites associated with exchange reactions was equal to 501. When I merged the the extracellular metabolites with the provided metabolite annotation file, only 137 were annotated. I have added all the results in the profile.
***************


## Files Added:
The models are in the models folder, which can be loaded by *load_models.m*.\

*parse_metabolite.py* parses the bigg metabolite database and extracts the metabolite id, name, Kegg and Chebi id only.\  

*exchange_metabolite.m* extracts metabolites associated with exchange reactions from all models provided and saves them in csv files.\   
*obtain_biggId.py* matches the exchange_metabolites names with the annotation file.\ 

 
