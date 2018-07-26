%Combined Model of microbiome, hepatic cell and small intestine epithelial cell(sIEC) mice. 
%References: 
% Liver metabolic reconstruction file: 
%               https://github.com/csbl/ratcon1
%               https://www.nature.com/articles/ncomms14250
% sIEC metabolic reconstruction file: 
%               https://www.nature.com/articles/s41598-017-07350-1
% microbiome models: 
%               4 models from SteadyCom exercise

% Read Models
livM = readCbModel(strcat('./models/','iRno.xml')); % liver model

sieC = load(strcat('./models/','mmu_ENT717.mat')); %iSEC model
sieC = sieC.model

% Microbiome models
cl =load(strcat('./models/','colstridium_ljungdahlii.mat'))
cl=cl.iHN637;
% convert the compartment format from e.g., '_c' to '[c]'
cl.mets = regexprep(cl.mets, '_([^_]+)$', '\[$1\]');
% make all empty cells in cell arrays to be empty string
fieldToBeCellStr = {'metFormulas'; 'genes'; 'grRules'; 'metNames'; 'rxnNames'; 'subSystems'};
for j = 1:numel(fieldToBeCellStr)
    cl.(fieldToBeCellStr{j})(cellfun(@isempty, cl.(fieldToBeCellStr{j}))) = {''};
end

bs = load(strcat('./models/','bsubtilis.mat'));
bs=bs.iYO844
% convert the compartment format from e.g., '_c' to '[c]'
bs.mets = regexprep(bs.mets, '_([^_]+)$', '\[$1\]');
% make all empty cells in cell arrays to be empty string
fieldToBeCellStr = {'metFormulas'; 'genes'; 'grRules'; 'metNames'; 'rxnNames'; 'subSystems'};
for j = 1:numel(fieldToBeCellStr)
    bs.(fieldToBeCellStr{j})(cellfun(@isempty, bs.(fieldToBeCellStr{j}))) = {''};
end


ll = load(strcat('./models/','lactococcus_lactis.mat'));
ll=ll.iNF517

% convert the compartment format from e.g., '_c' to '[c]'
ll.mets = regexprep(ll.mets, '_([^_]+)$', '\[$1\]');
% make all empty cells in cell arrays to be empty string
fieldToBeCellStr = {'metFormulas'; 'genes'; 'grRules'; 'metNames'; 'rxnNames'; 'subSystems'};
for j = 1:numel(fieldToBeCellStr)
    ll.(fieldToBeCellStr{j})(cellfun(@isempty, ll.(fieldToBeCellStr{j}))) = {''};
end


ec=load(strcat('./models/','iJO1366.mat'));
ec=ec.iJO1366
% convert the compartment format from e.g., '_c' to '[c]'
ec.mets = regexprep(ec.mets, '_([^_]+)$', '\[$1\]');
% make all empty cells in cell arrays to be empty string
fieldToBeCellStr = {'metFormulas'; 'genes'; 'grRules'; 'metNames'; 'rxnNames'; 'subSystems'};
for j = 1:numel(fieldToBeCellStr)
    ec.(fieldToBeCellStr{j})(cellfun(@isempty, ec.(fieldToBeCellStr{j}))) = {''};
end

