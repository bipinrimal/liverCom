models=modelCell(1:4,:)
modelNumber=size(models,1)

modelHosts=modelCell(5:6,:)
modelHostNum=size(modelHosts,1)


% Creating Name Tags
nameTagsModels={}
if isempty(nameTagsModels)
    for i = 1:modelNumber
        nameTagsModels{i,1}=strcat('model',num2str(i),'_');
    end
end

nameTagHost={}
if isempty(nameTagHost)
    for i=1:modelHostNum
        nameTagHost{i,1}=strcat('host',num2str(i),'_');
        
    end
end

for i = 1:modelNumber
    model=models{i,1}
    metIndices =~cellfun(@isempty,regexp(model.mets,'_e$'));
    model.mets(metIndices)=strrep(model.mets(metIndices),'_e','[e]')
    models{i,1}=model;
end

if ~isempty(modelHost)
    for i=1:modelHostNum
        modelHost=modelHosts{i,1}
        metIndices=~cellfun(@isempty,regexp(modelHost.mets,'_e$'));
        modelHost.mets(metIndices) = strrep(modelHost.mets(metIndices), '_e', '[e]');
        modelHosts{i,1}=modelHost
    end
end

eTag='u'
exTag='e'

%find fields in models (mets, metabolite Names)
presentinallModels=fieldnames(models{1})
missingFields={}
for i=2:modelNumber
    cfields=fieldnames(models{i})
    missingFields=union(missingFields,setxor(cfields,presentinallModels));
    missingFields
    presentinallModels=intersect(presentinallModels,cfields);
    size(presentinallModels)
end

models=restrictModelsToFields(models,presentinallModels)

modelStorage =cell(modelNumber,1)
for i=1:modelNumber
    model=models{i,1}
    exmod=model.rxns(strmatch('EX',model.rxns));
    model=removeRxns(model,exmod);
    modelStorage{i,1}=model;
end

