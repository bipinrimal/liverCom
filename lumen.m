nSp = 6;
spName = {'1','2','3','4','5','6'};
[rxnEx, metEx] = deal(cell(nSp, 1)); % matrix for exchange reactions
S = []; %stoichiometric matrix
metSps = []; %metabolites for species
rxnSps = []; % reactions
metCom2Sp = []; %common metabolites in species
rxnCom2Sp = []; %common reactions in species
[mets, rxns] = deal([]);

for j = 1:nSp
    rxnEx{j} = sum(modelCell{j}.S ~= 0, 1) == 1; % finds the exchange reactions, rows where exchange reactions are
    [metEx{j}, ~] = find(modelCell{j}.S(:, rxnEx{j})); % find the metabolite names corresponding the rxnEX
    [row, col] = size(S); %creates row and column of the size of S
    [rowJ, colJ] = size(modelCell{j}.S(:, ~rxnEx{j})); %creates row and columns of the size of exchange reactions
    S = [S                 sparse(row, colJ);...
     sparse(rowJ, col) sparse(modelCell{j}.S(:, ~rxnEx{j}))]; %filling up?
    metSps = [metSps; [j * ones(rowJ, 1), (1:rowJ)']]; % indices for metabolites
    rxnSps = [rxnSps; [j * ones(colJ, 1), find(~rxnEx{j}(:))]];
    mets = [mets; strcat(modelCell{j}.mets, '_', spName{j})];
    rxns = [rxns; strcat(modelCell{j}.rxns(~rxnEx{j}), '_', spName{j})];
end
% add the lumen compartment
% get a unique list of lumen metabolites
metU = {}; % empty lumen
for j = 1:nSp
    metU = union(metU, modelCell{j}.mets(metEx{j})); %.metCom (union of all metabolites not good..need to be unique)
end
S = [S; sparse(numel(metU), size(S, 2))]; %adding lumen 
metSps = [metSps; [(nSp + 1) * ones(numel(metU), 1), (1:numel(metU))']]; % indices
metU(cellfun(@isempty, metU)) = [];
mets = [mets; strcat(metU, '_u')];

nRxnU = 0;
for j = 1:nSp
    [yn, id] = ismember(modelCell{j}.mets, metU); %.metCom Finds if metU elements are in modelCell{j}.mets
    
    for k = 1:numel(modelCell{j}.mets) %for all metabolites in modelCell{j}
        if yn(k)
            nRxnU = nRxnU + 1;
            % row number for the met in cell j
            rK = find(metSps(:, 1) == j & metSps(:, 2) == k);
            % row number for the met in the lumen
            rU = find(metSps(:, 1) == nSp + 1 & metSps(:, 2) == id(k));
            S = [S, sparse([rK; rU], [1; 1], [-1; 1], size(S, 1), 1)];
            rxns = [rxns; ['EX_' metU{id(k)} '_' spName{j} '_u']];
            rxnSps = [rxnSps; (nSp + 1), nRxnU];
        end
    end
end

%Exchange between liver and iSEC:

