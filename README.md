# liverCom
microbiome-liver model

The model incorporates the microbiome community and its interaction with host cell; mainly hepatocyte through the intestinal epithelial cell. In preliminary work, there were 4 bacteria used Bacillus subtilis, Colstridium ljundahlii, Lactococcus lactii, and Escherichia coli. The host models used for the integrated model involved intestinal epithelial cell model of mice and thus the liver model was also of mice to keep it updated. 

1. The first script *load_models.m* loads the models (both sbml and mat formats). The metabolites compartment names in each model are converted from ‘_c’ to ‘[c]’ standard. All empty cells in cell arrays are converted to empty string. 

The metabolites in intestinal epithelial cell model has a single ‘_’ compared to ‘__’ in bacterial models. This was rectified. 

The liver model has non-standard naming system for the metabolites.

2. The second script *liver_lumen_exchanges.m* tries to address the non-conventional naming of metabolites. Manually, going through the exchange reactions of both intestinal epithelial cell and the liver cell, the same metabolites (gleaned from the metabolite names) are converted to common name. The file is *common_iseCvslivM.xls* Here, the liver metabolite names are converted to match intestinal epithelial cell. 

3. The third script *lumen.m* is the script which generates the combined model. There are three parts to the script:
a)	First script generates the giant matrix with all the internal reactions.
For each model, the exchange reactions are identified using the stoichiometric coefficients. 

For liver model, we are only concerned on what is coming out of the model into the lumen that is relevant to us. Those are bile acids. Going through the model again manually, all those reactions involving bile acids are selected and all the other reactions are considered internal reactions. 

b)	The second script introduces a common lumen compartment for the exchange of metabolites between the models. 
First, the union of the exchange metabolites for the models (except the liver) is calculated. Union is not ideal for this type of compartment building as each should be unique. However, manually going through the metU metabolites, their uniqueness was determined. 

A list of rows for the unique metabolites was added to the giant matrix and based on whether these metabolites are in the microbiome or intestinal model, pertinent reactions are added to the matrix. Each reaction is named accordingly. 

c)	The third script adds an exchange compartment for the transfer of metabolites from intestinal epithelial cell to the liver cell. Instead of another compartment, adding an exchange reaction from the intestinal epithelial cell to the liver cell is desired. However, for the ease of making it, another compartment was added. 

In this compartment, the metabolites coming out of the intestinal epithelial cell and also present in the liver cell are exchanged. Hence, an intersection of the common metabolites gives a list of metabolites common to both. The reversibility can be addressed later to allow only one directional flow of metabolites. 

d)	The final script is adding the bile acid exchange compartment. 
The bile acids are deposited into the lumen. Any bacteria that has the potential to take up the bile acids will import the bile acid produced by the liver. Thus, the script goes through the metabolites of each bacterial cell and checks whether the metabolite matches the bacterial metabolites and adds the reaction to the matrix. 

The final script merges the giant S matrix, the reactions, the metabolites, the ub, lb, b and c calculated in the script. 


Problems: 
1.	The bacteria involved in the model do not have the bile acids produced by the liver which means it is just a giant model without the liver model. I am using the 9-bacteria model produced in the paper to see if there are bile acids in there. 

2.	I do not have a biomass objective for the model. 

3.	As you guess, optimize Cb Model leads to unfeasible result.


Debolina suggested using a single identifier for all metabolites. For example; Kegg identifier. That did not work. A lot of the metabolites in liver models (even for bile acids were ambiguous). So that did not work. 

