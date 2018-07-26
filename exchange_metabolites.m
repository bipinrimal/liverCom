% find exchange reactions for all
modelCell=struct('bs',bs,'cl',cl,'ec',ec,'ll',ll,'livM',livM,'sieC',sieC)
spNames = fieldnames(modelCell)

%confvert to cell array for indexing
modelCell=struct2cell(modelCell) 
nSp = numel(modelCell);


% Needed
rxnEx=cell(nSp,1)
ex = cell(nSp,1)
nameList=cell(nSp,1)
%loop through the structures to get exchange Metabolites
for j=1:nSp
    % Find the columns in stoichiometric matrix with only one non-zero
    % value (+1,-1) indicating exchanging reactions. 
    rxnEx{j}=sum(modelCell{j}.S~=0,1)==1 % ~=0 (not equal to zero); 1 means column, ==1 means only one value)
    ex{j}=any(modelCell{j}.S(:,rxnEx{j}),2)
    nameList{j}=modelCell{j}.mets(ex{j})
    %filename
    filename=strcat(spNames{j},'_exMet.csv')
    %save file
    fid=fopen(filename,'wt');
    if fid>0
        for k=1:size(nameList{j},1)
            fprintf(fid,'%s,\n',nameList{j}{k,:});
        end
        fclose(fid);
    end
end
    
