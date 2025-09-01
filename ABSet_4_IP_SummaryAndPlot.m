
%% confusion matrix
OutPath=['E:\PFC and MTL\Results\Review\decoding\ABSet_4\IP_Bar\'];
GroupName = {'lPFC','mPFC','HPC'};
CellScreen = {'','_3','_correct','_correct_3','_IP_4'};
nCellScreen = 5;
%% all possible pairs
Decoder=1;
for nregion = 1:3
    result_path = [OutPath,'\',GroupName{nregion},'_all_AllComb_900\',CellScreen{nCellScreen},'\'];
    load([result_path,'\',GroupName{nregion},'_all.mat'],'PredictFit')
    PredictFit_All=cell(2,1);
    for ncomb=1:size(PredictFit,1)
        for nrep=1:size(PredictFit,2)
            for npre=1:2
                PredictFit_All{npre}=cat(1,PredictFit_All{npre},PredictFit{ncomb,nrep,npre});
            end
        end
    end
    clear PredictFit
    subplot(1,3,nregion)
    decodedLabels=PredictFit_All{1};
    trueLabels=PredictFit_All{2}(:,Decoder(1));
    tool_ConfusionMatrix_v2(decodedLabels,trueLabels)
    % customMap=tool_WinterColorMap(0.25);
     customMap=tool_BiColorMap(0.25);
    colormap(customMap)
    title(['Co-location decoding',GroupName{nregion}])
end
pause
FilePath=[OutPath,'/ABset_4.jpg'];
saveas(gcf,FilePath)
close all

