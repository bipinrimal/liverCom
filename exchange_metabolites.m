% find exchange reactions for all
modelCell=struct('bs',bs,'cl',cl,'ec',ec,'ll',ll,'sieC',sieC,'livM',livM)
spNames = fieldnames(modelCell)

%confvert to cell array for indexing
modelCell=struct2cell(modelCell) 
nSp = numel(modelCell);


% Needed
rxnEx=cell(nSp,1)
ex = cell(nSp,1)
commonNames = cell(nSp,1)
nameList=cell(nSp,1)
pos=cell(nSp,1)
name_index=cell(nSp,1)
%loop through the structures to get exchange Metabolites
for j=1:nSp
    % Find the columns in stoichiometric matrix with only one non-zero
    % value (+1,-1) indicating exchanging reactions. 
    rxnEx{j}=sum(modelCell{j}.S~=0,1)==1 % ~=0 (not equal to zero); 1 means column, ==1 means only one value)
    ex{j}=any(modelCell{j}.S(:,rxnEx{j}),2)
    nameList{j}=modelCell{j}.mets(ex{j})
    pos{j}=string(find(ex{j}))
    commonNames{j}=modelCell{j}.metNames(ex{j})
    met_formulas{j}=modelCell{j}.metFormulas(ex{j})
    name_index{j} = [commonNames{j},pos{j},nameList{j}]
    %filename
%     filename=strcat(spNames{j},'_exMet.csv')
%     comName=strcat(spNames{j},'_exCommonNames.csv')
    posName=strcat(spNames{j},'_exIndices.csv')
    %save file
%     fid=fopen(filename,'wt');
%     fCid=fopen(comName,'wt');
%     fkid=fopen(posName,'wt');
%     
%     cell2csv(filename,nameList{j},'\t')
%     cell2csv(comName,commonNames{j},'\t')
    cell2csv(posName,name_index{j},'\t')
%     if fid>0
%         for k=1:size(nameList{j},1)
%             fprintf(fid,'%s,\n',nameList{j}{k,:});
%         end
%         
%         fclose(fid);
%     end
%     if fCid>0
%         for k=1:size(commonNames{j},1)
%             fprintf(fCid,'%s,\n',commonNames{j}{k,:});
%         end
%         
%         fclose(fCid);
%     end
%     
%     if fkid>0
%         for k=1:size(name_index{j},1)
%             fprintf(fkid,'%s,\n',name_index{j}{k,:});
%         end
%         
%         fclose(fkid);
%     end
%     
end
    
%% Converting to standard
%Read the kegg api database file
%kegg_db=readtable('kegg_api.xls','ReadRowNames',true);
% for each model
%common_metabolites={}
% for j=1:nSp
%     %loop through the exchange metabolites
%     for k=metEx{j}
%         %search for the metabolite name in the kegg_db
%         pattern=modelCell{j}.metNames{k};     
%         search=contains(kegg_db.name,pattern,'IgnoreCase',true);
%  
%         %take out the compartment information
%         %compartment=regexp(modelCell{j}.mets{k}, '{3}$','match')
%         metabolite=modelCell{j}.mets{k}
%         compartment=metabolite(length(metabolite)-2:end)
%         
%         %if there is only one match
%         if sum(search)==1
%             met_id=kegg_db.kegg_id(search);
%             met_id=erase(met_id,'cpd:');
%             kegg_id=strcat(met_id,compartment);
%             modelCell{j}.keggs{k}=strcat(met_id,compartment);
%             
%         if sum(search)==0
%             modelCell{j}.keggs{k}=modelCell{j}.mets{k};
%             
%         %if there is more than one match, check if the search in kegg database has ";" at the end of pattern.    
%         if sum(search)>1
%             for l=1:sum(search)
%                 multiple_names=kegg_db.kegg_name(search);
%                 multiple_ids =kegg_db.kegg_id(search);
%                 index_of_colon=regexp(multiple_names{l},pattern)+length(pattern);
%                 if multiple_names{l}{index_of_colon}==';'
%                     met_id=multiple_ids{l};
%                     met_id=erase(met_id,'cpd:');
%                     modelCell{j}.keggs{k}=strcat(met_id,compartment);
%                 else
%                     continue
%                 end
%             end
%         end
%         end
%         end
%     end
% end


%% Converting to standard
%Read the kegg api database file
kegg_db=readtable('kegg_api.xls');
% for each model
for j=1:nSp
    [yn,id]=ismember(modelCell{j}.mets,modelCell{j}.mets(metEx{j})); %the indices for the exchange mets
    %loop through the exchange metabolites
    for k = 1:numel(modelCell{j}.mets)
        if yn(k)
            %search for the metabolite name in the kegg_db
            pattern=modelCell{j}.metNames{k};     
            search=contains(kegg_db.name,pattern,'IgnoreCase',true);

            %take out the compartment information
            compartment=regexp(modelCell{j}.mets{k}, '[.\]','match');

            %if there is only one match
            if sum(search)==1
                met_id=kegg_db.kegg_id(search);
                met_id=erase(met_id,'cpd:');
                kegg_id=strcat(met_id,compartment);
                modelCell{j}.keggs{k}=strcat(met_id,compartment);

            if sum(search)==0
                modelCell{j}.keggs{k}=modelCell{j}.mets{k};

            %if there is more than one match, check if the search in kegg database has ";" at the end of pattern.    
            if sum(search)>1
                for l=1:sum(search)
                    multiple_names=kegg_db.name(search);
                    multiple_ids =kegg_db.kegg_id(search);
                    index_of_pattern=regexp(multiple_names{l},pattern);
                    if multiple_names{l}{index_of_pattern+length(pattern)}==';'
                        if multiple_names{l}{index_of_pattern - 1}==' '
                            met_id=multiple_ids{l};
                            met_id=erase(met_id,'cpd:');
                            modelCell{j}.keggs{k}=strcat(met_id,compartment);
                        else
                            modelCell{j}.keggs{k}=modelCell{j}.mets{k};
                        end
                    else
                        modelCell{j}.keggs{k}=modelCell{j}.mets{k};
                    end
                end
            end
            end
            end
        end
    end
end


    
