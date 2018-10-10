zh54%Read thze models
cl =load('colstridium_ljungdahlii.mat')
cl=cl.iHN637;
% convert the compartment format from e.g., '_c' to '[c]'
cl.mets = regexprep(cl.mets, '_([^_]+)$', '\[$1\]');
% make all empty cells in cell arrays to be empty string
fieldToBeCellStr = {'metFormulas'; 'genes'; 'grRules'; 'metNames'; 'rxnNames'; 'subSystems'};
for j = 1:numel(fieldToBeCellStr)
    cl.(fieldToBeCellStr{j})(cellfun(@isempty, cl.(fieldToBeCellStr{j}))) = {''};
end

bs = load('bsubtilis.mat');
bs=bs.iYO844
% convert the compartment format from e.g., '_c' to '[c]'
bs.mets = regexprep(bs.mets, '_([^_]+)$', '\[$1\]');
% make all empty cells in cell arrays to be empty string
fieldToBeCellStr = {'metFormulas'; 'genes'; 'grRules'; 'metNames'; 'rxnNames'; 'subSystems'};
for j = 1:numel(fieldToBeCellStr)
    bs.(fieldToBeCellStr{j})(cellfun(@isempty, bs.(fieldToBeCellStr{j}))) = {''};
end


ll = load('lactococcus_lactis.mat');
ll=ll.iNF517

% convert the compartment format from e.g., '_c' to '[c]'
ll.mets = regexprep(ll.mets, '_([^_]+)$', '\[$1\]');
% make all empty cells in cell arrays to be empty string
fieldToBeCellStr = {'metFormulas'; 'genes'; 'grRules'; 'metNames'; 'rxnNames'; 'subSystems'};
for j = 1:numel(fieldToBeCellStr)
    ll.(fieldToBeCellStr{j})(cellfun(@isempty, ll.(fieldToBeCellStr{j}))) = {''};
end


ec=load('iJO1366.mat');
ec=ec.iJO1366
% convert the compartment format from e.g., '_c' to '[c]'
ec.mets = regexprep(ec.mets, '_([^_]+)$', '\[$1\]');
% make all empty cells in cell arrays to be empty string
fieldToBeCellStr = {'metFormulas'; 'genes'; 'grRules'; 'metNames'; 'rxnNames'; 'subSystems'};
for j = 1:numel(fieldToBeCellStr)
    ec.(fieldToBeCellStr{j})(cellfun(@isempty, ec.(fieldToBeCellStr{j}))) = {''};
end


%Sanity check of models



%polishing the model




nameTagsModel = {'bs';'col';'ec';'ll'};
EcCom = createMultipleSpeciesModel({bs; cl; ec; ll}, nameTagsModel);
EcCom.csense = char('E' * ones(1,numel(EcCom.mets)));  % correct the csense


[EcCom.infoCom, EcCom.indCom] = getMultiSpeciesModelId(EcCom, nameTagsModel);
%disp(EcCom.infoCom);

%Biomass Reactions: 
%bs:{'BIOMASS_BS_10'}
%col: {'BIOMASS_Cl_DSM_WT_46p666M1'}
%ec: {'BIOMASS_Ec_iJO1366_WT_53p95M'}
%ll: {'BIOMASS_LLA'}

biomassRxns = [{'BIOMASS_BS_10'} {'BIOMASS_Cl_DSM_WT_46p666M1'} 
    {'BIOMASS_Ec_iJO1366_WT_53p95M'} {'BIOMASS_LLA'}];

for j = 1:numel(biomassRxns)
    rxnBiomass = strcat(nameTagsModel,biomassRxns{j});
    rxnBiomassId=findRxnIDs(EcCom,rxnBiomass);
    EcCom.infoCom.spBm=rxnBiomass;
    EcCom.indCom.spBm=rxnBiomassId;
end


models=[{bs} {cl} {ec} {ll}]
for j=1:numel(models)
    metabolites=models{j}.mets
    metEx=strcmp(getCompartment(metabolites),'e');
    
    rxnExAll=find(sum(models{j}.S~=0,1)==1);
    [rxnEx,~]=find(models{j}.S(metEx,rxnExAll)');
    rxnEx=rxnExAll(rxnEx);
    lbEx=models{j}.lb(rxnEx);
    
    [yn, id] = ismember(strrep(models{j}.mets(metEx), '[e]', '[u]'), ...
        EcCom.infoCom.Mcom);  % map the metabolite name
    %assert(all(yn));  % must be a 1-to-1 mapping
    EcCom.lb(EcCom.indCom.EXcom(:,1)) = lbEx(id);  % assign community uptake bounds
    EcCom.ub(EcCom.indCom.EXcom(:,1)) = 1e5;
    f = EcCom.indCom.EXsp(:, j) ~= 0;
    EcCom.lb(EcCom.indCom.EXsp(f, j))=lbEx(id);  % assign organism-specific uptake bounds
end
