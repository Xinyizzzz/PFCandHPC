%% plot for each set
clear all
clc
close all
%% confusion matrix
OutPath=['E:\PFC and MTL\Results\Review\decoding\ItemIdentity\IP_Bar\'];
GroupName = {'lPFC','mPFC','HPC'};
CellScreen = {'','_3','_correct','_correct_3','_IP_4'};
nCellScreen = 5;
%% all possible pairs
Decoder=1;
kk=0;
for nregion = 1:3
    result_path = [OutPath,'\',GroupName{nregion},'_all_5\',CellScreen{nCellScreen},'\'];
    load([result_path,'\',GroupName{nregion},'_all.mat'],'PredictFit')
    PredictFit_All=cell(4,2);
    for npair=1:4
        for nrep=1:size(PredictFit,2)
            for npre=1:2
                PredictFit_All{npair,npre}=cat(1,PredictFit_All{npair,npre},PredictFit{npair,nrep,npre});
            end
        end
    end
    for npair=1:4
        kk=kk+1;
        subplot(3,4,kk)
        decodedLabels=PredictFit_All{npair,1};
        trueLabels=PredictFit_All{npair,2}(:,Decoder(1));
        tool_ConfusionMatrix_v2(decodedLabels,trueLabels)
        title(['item identity set',num2str(npair),' ',GroupName{nregion}])
    end
end
pause
FilePath=[OutPath,'/ItemIdentity_each.jpg'];
saveas(gcf,FilePath)
close all


%% average across all set
clear all
clc
close all
%% confusion matrix
OutPath=['E:\PFC and MTL\Results\Review\decoding\ItemIdentity\IP_Bar\'];
GroupName = {'lPFC','mPFC','HPC'};
CellScreen = {'','_3','_correct','_correct_3','_IP_4'};
nCellScreen = 5;
%% all possible pairs
Decoder=1;
kk=0;
figure()
set(gcf,'OuterPosition',[992 248 1236 1075])
for nregion = 1:3
    result_path = [OutPath,'\',GroupName{nregion},'_all_5\',CellScreen{nCellScreen},'\'];
    load([result_path,'\',GroupName{nregion},'_all.mat'],'PredictFit')
    PredictFit_All=cell(2,1);
    for npair=1:4
        for nrep=1:size(PredictFit,2)
            for npre=1:2
                PredictFit_All{npre}=cat(1,PredictFit_All{npre},PredictFit{npair,nrep,npre});
            end
        end
    end
    kk=kk+1;
    subplot(1,3,kk)
    decodedLabels=PredictFit_All{1};
    trueLabels=PredictFit_All{2}(:,Decoder(1));
    tool_ConfusionMatrix_v2(decodedLabels,trueLabels)
    customMap=tool_BiColorMap(0.5);
    colormap(customMap)
    colorbar
    title(['Item identity ',GroupName{nregion}],'FontWeight','normal')
                set(gca,'FontSize',18)
end
pause
FilePath=[OutPath,'/ItemIdentity_All.jpg'];
saveas(gcf,FilePath)
close all

