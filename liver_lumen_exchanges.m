modelCell=struct('bs',bs,'cl',cl,'ec',ec,'ll',ll,'sieC',sieC,'livM',livM)
spNames = fieldnames(modelCell)

%convert to cell array for indexing
modelCell=struct2cell(modelCell) 
nSp = numel(modelCell);
%% import the exchange_indices table
exchange_indices=readtable('common_iseCvs livM.xls');

% loop through the table
for k = 1:size(exchange_indices, 1)
    liver_index=table2array(exchange_indices(k,2));
    iseC_index=table2array(exchange_indices(k,4));
    %liver_metName=modelCell{6}.mets(liver_index);
    iseC_metName=modelCell{5}.mets(iseC_index);
    %replace compartment
    %iseC_metName=regexprep(iseC_metName,'\[.\]','\[s\]')
    modelCell{6}.mets(liver_index)=iseC_metName 
end 
