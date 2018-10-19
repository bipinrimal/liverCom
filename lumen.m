nSp = 6;
spName = {'1','2','3','4','5','6'};
[rxnEx, metEx] = deal(cell(nSp, 1)); % matrix for exchange reactions
S = []; %stoichiometric matrix
%rev = [];
lb = [];
ub = [];
c = [];
b = [];
%map rxn
metSps = []; %metabolites for species
rxnSps = []; % reactions
metCom2Sp = []; %common metabolites in species
rxnCom2Sp = []; %common reactions in species
[mets, rxns] = deal([]);
rxnEx2met = cell(nSp, 1);

%% A giant matrix with all the internal reactions. You create this matrix and forget about it
for j = 1:nSp
    
%     lbJ = modelCell{j}.lb;
%     ubJ = modelCell{j}.ub;
%     cJ = modelCell{j}.c;
    
    %## this is sufficient for finding all exchange reactions, unless
    %## you want to exclude some?
    %## I reverted your change for bile acids. See the code after this loop
    rxnEx{j} = (sum(modelCell{j}.S ~= 0, 1) == 1)';
    %& (sum(modelCell{j}.S(~rxnEx{j},:) ~= 0) == 0)';
    
   
    
    %metLmets = find(modelCell{6}.S(:, bile_exRxns))
    % finds the exchange reactions,
    %rows where exchange reactions are 
    %(I think just the 4 microbes and the intestinal epithelial cell, and
    %liver goes different)edit:2 does not matter because this is for
    %internal reactions
    [metEx{j}, ~] = find(modelCell{j}.S(:, rxnEx{j})); % find the metabolite names corresponding the rxnEX
    [row, col] = size(S); %creates row and column of the size of S
    [rowJ, colJ] = size(modelCell{j}.S(:, ~rxnEx{j})); %creates row and columns of the size of internal reactions
    S = [S                 sparse(row, colJ);...
     sparse(rowJ, col) sparse(modelCell{j}.S(:, ~rxnEx{j}))]; %filling up internal reactions of everyone
    metSps = [metSps; [j * ones(rowJ, 1), (1:rowJ)']]; % indices for metabolites, count 1 to number or rowJ, and create indices with species identifier j
    rxnSps = [rxnSps; [j * ones(colJ, 1), find(~rxnEx{j}(:))]]; 
    mets = [mets; strcat(modelCell{j}.mets, '_', spName{j})];
    rxns = [rxns; strcat(modelCell{j}.rxns(~rxnEx{j}), '_', spName{j})];
    
    for k = 1:colJ
        if rxnEx{j}(k)
            metJK = find(modelCell{j}.S(:, k), 1);
            rxnEx2met{j}(metJK,:) = [col + k, k];
            %set positive flux of exchange reaction as uptake
            %(convention)
            if modelCell{j}.S(metJK, k) > 0
                s = modelCell{j}.S(metJK,k);
                modelCell{j}.S(metJK,k) = -1;
                [ubJ(k), lbJ(k), cJ(k)] = deal(-s * lbJ(k), -s * ubJ(k), cJ(k));
            end
        end
    end
    %rev = [rev; modelCell{j}.rev];
    lbJ = modelCell{j}.lb(~rxnEx{j});
    ubJ = modelCell{j}.ub(~rxnEx{j});
    cJ = modelCell{j}.c(~rxnEx{j});
    
    lb = [lb; lbJ];
    ub = [ub; ubJ];
    c = [c; cJ];
    b = [b; modelCell{j}.b];
    
    total_internal=length(rxns)
    total_lb = length(lb)

end

% The bile acid exchange reactions in the model (manual curation)
%## do not mess the bile acid exchange directly with all exchange
%## reactions. You still need sink/source for other extracellular
%## metabolites otherwise the model cannot grow. So treat bile acids as
%## some special objects here
%## better use the reaction ids than the indices, which may change if
%## you modify the model beforehand
tmp = {'RCR30067';'RCR30473';'RCR30474';'RCR90070';'RCR90071';'RCR90072';'RCR90073';'RCR90074';'RCR90075';'RCR90077';'RCR90080';'RCR90082';'RCR90083';'RCR90084';'RCR90085';'RCR90086';'RCR90087';'RCR90088';'RCR90089';'RCR90090';'RCR90091';'RCR90092';'RCR90093';'RCR90094';'RCR90179'};
rxnBile = false(numel(modelCell{6}.rxns), 1);
rxnBile(findRxnIDs(modelCell{6}, tmp)) = true;
[metBileId, ~] = find(modelCell{6}.S(:, rxnBile));
metBile = modelCell{6}.mets(metBileId);
metBileLogic = false(numel(modelCell{6}.mets), 1);
metBileLogic(metBileId) = true;
%rxnEx{6} = [5963,6357, 6358,8111,8112,8113,8114,8115,8116,8118,8121,8123, ...
%           8124,8125,8126,8127,8128,8129,8130,8131,8132,8133,8134,8135,8205];

% All non-exchange reactions are wrong here. For liver, not all exchange
% reactions are in the lumen


%% The Lumen
% add the lumen compartment
% get a unique list of lumen metabolites

metU = {}; % empty lumen

for j = 1:5
    metU = union(metU, modelCell{j}.mets(metEx{j})); %.metCom (union of all metabolites not good..need to be unique)
end

%## add bile acids into metU
metU = union(metU, metBile);

S = [S; sparse(numel(metU), size(S, 2))]; %adding lumen compartment as another matrix of rows below the stoichiometric matrix
metSps = [metSps; [(nSp + 1) * ones(numel(metU), 1), (1:numel(metU))']]; % indices for the added common lumen metabolites
metU(cellfun(@isempty, metU)) = [];
mets = [mets; strcat(metU, '_u')];
b = [b; zeros(numel(metU),1)];
nRxnU = 0;
for j = 1:5
    [yn, id] = ismember(modelCell{j}.mets, metU); %.metCom Finds if metU elements are in modelCell{j}.mets
    
    for k = 1:numel(modelCell{j}.mets) %for all metabolites in modelCell{j}
        if yn(k) % if logical value for metabolite number is true, yn(k) is logical value for number k obtained from above, if its true, it runs
            nRxnU = nRxnU + 1; %simple counter
            % row number for the met in cell j
            rK = find(metSps(:, 1) == j & metSps(:, 2) == k);
            % row number for the met in the lumen
            rU = find(metSps(:, 1) == nSp + 1 & metSps(:, 2) == id(k));
            S = [S, sparse([rK; rU], [1; 1], [-1; 1], size(S, 1), 1)];
            rxns = [rxns; ['EX_' metU{id(k)} '_' spName{j} '_u']];
            rxnSps = [rxnSps; (nSp + 1), nRxnU];
            ub=[ub;10000 * ones];
            %## I would say in general you should leave all exchange
            %## bidirection unless you want to limit the model to have only a
            %## certain number of exchanges possible
            lb = [lb; -10000];
            %lb = [lb; zeros];
            c = [c; zeros];
        end
    end
end


% ub = [ub;10000 * ones(numel(metU),1)];
% lb = [lb; zeros(numel(metU),1)];
% c = [c;  zeros(numel(metU),1)];
% %rev = [rev;  zeros(2*numel(metU),1)];
lumen_rxns=length(rxns)
lumen_lb = length(lb)


%% Exchange compartment (blood stream)

%## It is better to do the same here as in the lumen. And if needed, limit
%## the exchange between the blood stream and the organ (e.g., no bile acids
%## into blood stream). It will be easier if you want to add new organs in the
%## future
%## I also suggest using 'b' for blood instead of 'i', but completely up to you
metI={};
metI = union(modelCell{6}.mets(metEx{6}), modelCell{5}.mets(metEx{5}));  %metabolites that are exchange reactions in Intestinal Epithelial Cell
metI(cellfun(@isempty,metI))=[];
S = [S;sparse(numel(metI),size(S,2))];%adding exchange compartment between IseC and liver cell
metSps=[metSps; [(nSp+2)*ones(numel(metI),1),(1:numel(metI))']];
mets = [mets; strcat(metI, '_i')];
b = [b; zeros(numel(metI),1)];
nRxnI=0;
% metabolites to be excluded from exchanges with blood
excludeFromBlood = cell(6, 1);
% nothing for IEC
excludeFromBlood{5} = false(numel(modelCell{5}.mets), 1);
% do not add exchange into blood stream for bile acids (add others if you know)
excludeFromBlood{6} = metBileLogic;
for j=5:6
    [yn,id]=ismember(modelCell{j}.mets, metI);
    for k=1:numel(modelCell{j}.mets)
        if yn(k) && ~excludeFromBlood{j}(k)
            nRxnI=nRxnI+1;
            rk=find(metSps(:,1)==j & metSps(:,2)==k);
            rI=find(metSps(:,1)==nSp + 2 & metSps(:,2)==id(k));
            S=[S,sparse([rk;rI],[1;1],[-1,1],size(S,1),1)];
            rxns=[rxns;['EX_' metI{id(k)} '_' spName{j} '_i']];
            rxnSps=[rxnSps;(nSp+2),nRxnI];
            ub=[ub; 10000];
            lb = [lb; -10000];
            c = [c; 0];
        end
    end
end

% ub = [ub;10000 * ones(numel(metI),1)];
% lb = [lb; zeros(numel(metI),1)];
% c = [c;  zeros(numel(metI),1)];
% exchange_rxns = length(rxns)
% exchange_lb = length(lb)
%rev = [rev;  zeros(2*numel(metI),1)];


%% Bile acid exchange

% the bile acids are deposited into the lumen. 
% The bile acid exchange reactions in the model (manual curation)

% Find the metabolites associated with the bile acid exchange reactions
% ..i.e. exchange bile acids metabolites

%metL=modelCell{6}.mets(metBileId);

%## I think a better way here is just to include the bile acids as
%## metabolites in the lumen and then just add exchange reactions from liver to lumen. 
%## No need to add exchange reactions for each bacteria and no need to add
%## an additional compartment

% index of metBile in metU
[yn, id] = ismember(metBile, metU);
for j = 1:numel(metBile)
    % find where bile acid is in the lumen and the row
    rk=find(metSps(:,1) == nSp + 1 & metSps(:,2) == id(j));
    rL=find(metSps(:,1) == 6 & metSps(:,2) == metBileId(j));
    % add the stoichiometric so that the exchange can happen
    S = [S,sparse([rk;rL],[1;1],[1, -1],size(S,1),1)];
end
% name the reactions
rxns = [rxns; strcat('EX_', metBile, '_', spName{6}, '_l')];
rxnSps  =[rxnSps; (nSp+3) * ones(numel(metBile), 1), (1:numel(metBile))'];
ub = [ub; 10000 * ones(numel(metBile), 1)];
% these will be irreversible exchanges because only from liver to lumen
lb = [lb; zeros(numel(metBile), 1)];
c = [c; zeros(numel(metBile), 1)];

% %create a exchange platform for bile aicds below the giant matrix
% S=[S;sparse(numel(metL),size(S,2))];
% b = [b; zeros(numel(metL),1)];
% 
% %the indices for the metabolites added to the end of the giant matrix
% metSps=[metSps;[(nSp+3)*ones(numel(metL),1),(1:numel(metL))']];
% %metL(cellfun(@isempty,metL))=[]
% mets = [mets; strcat(metL, '_l')];
% 
% 
% %then once they are in the lumen, if any of the bacteria (or iSEC) can take in the
% %bile acid, they should be taken in. Thus, running the loop for all except liver and Isec 
% nRxnL=0
% for j=1:4
%     % Are either of the metL metabolites in the bacteria and ISEC?
%     [yn,id]=ismember(modelCell{j}.mets, metL);
%     %for all metabolites
%     for k=1:numel(modelCell{6}.mets)
%         if yn(k) % if bile acid present in the bacteria
%             nRxnL=nRxnL+1; % count the reaction
%             % find where bile acid is in the model of bacteria and the row
%             rk=find(metSps(:,1)==j & metSps(:,2)==k);
%             % find where the bile acid is in the exchange compartment
%             rL=find(metSps(:,1)==nSp+3 & metSps(:,2)==id(k));
%             % add the stoichiometric so that the exchange can happen
%             S=[S,sparse([rk;rI],[1;1],[-1,1],size(S,1),1)];
%             % name the reaction
%             rxns=[rxns;['EX_',metL{id(k)} '_',spName{6} '_l']];
%             rxnSps=[rxnSps;(nSp+3),nRxnL];
%             ub=[ub;10000 * ones];
%             lb = [lb; zeros];
%             c = [c; zeros];
%         end
%     end
% end
% % bile_rxns=length(rxns)
% % bile_lb = length(lb)
% % ub = [ub; 10000 * ones(numel(metL),1)];
% % lb = [lb; zeros(numel(metL),1)];
% % c = [c;  zeros(numel(metL),1)];
% %rev = [rev;  zeros(2*numel(metL),1)];


%% system input/output reactions
%## You need to add them for the model to be really able to have any uptake
%## from/to the environment (here are the diet/excretion in the lumen (metU) 
%## and plasma metabolites in the blood stream (metI)

metUid = find(metSps(:, 1) == nSp + 1);
S = [S, sparse(metUid, 1:numel(metU), -1, size(S,1), numel(metU))];
rxns = [rxns; strcat('EX_', metU, '_u')];
rxnSps = [rxnSps; (nSp + 1) * ones(numel(metU), 1), nRxnU + (1:numel(metU))'];
ub=[ub; 10000 * ones(numel(metU), 1)];
lb = [lb; zeros(numel(metU), 1)];
c = [c; zeros(numel(metU), 1)];

%## define the nutrient available (from diet) in the lumen using the lb
%## for the exchange reactions above

metIid = find(metSps(:, 1) == nSp + 2);
S = [S, sparse(metIid, 1:numel(metI), -1, size(S,1), numel(metI))];
rxns = [rxns; strcat('EX_', metI, '_i')];
rxnSps = [rxnSps; (nSp + 2) * ones(numel(metI), 1), nRxnI + (1:numel(metI))'];
ub=[ub; 10000 * ones(numel(metI), 1)];
lb = [lb; zeros(numel(metI), 1)];
c = [c; zeros(numel(metI), 1)];

%## define the nutrient available in the blood stream (e.g., O2, glucose, lactate)
%## using the lb for the exchange reactions above

