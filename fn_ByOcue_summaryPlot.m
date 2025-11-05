% function fn_ByOcue_summaryPlot
GroupName = {'lPFC','mPFC','HPC'};
CellScreen = {'','_3','_correct','_correct_3'};
nCellScreen = 3;
%% all possible pairs
Pair={'-90','0';'-90','90';'0','-90';'0','90';'90','0';'90','-90'};
for nregion = 1:3
    for nfolder = 1:length(Pair(:,1))
        result_path = ['E:\PFC and MTL\Results\Review\decoding\ByOcue\Bar_OP\',GroupName{nregion},'_all\',Pair{nfolder,1},'Train',Pair{nfolder,2},'Test',CellScreen{nCellScreen},'\'];
        load([result_path,'\',GroupName{nregion},'_all.mat'],'performanceRER')
        Aved(nregion,1,nfolder)=performanceRER(1,1);
        Aved(nregion,2,nfolder)=performanceRER(4,1);
    end
end
Avem=mean(Aved,3);

% % function fn_ByOcue_summaryPlot
% GroupName = {'lPFC','mPFC','HPC'};
% CellScreen = {'','_3','_correct','_correct_3'};
% nCellScreen = 3;
% %% all possible pairs
% Pair={'-90','0';'-90','90';'0','-90';'0','90';'90','0';'90','-90'};
% for nregion = 1:3
%     for nfolder = 1:length(Pair(:,1))
%         result_path = ['E:\PFC and MTL\Results\decoding\ByOcue\Bar_OP\',GroupName{nregion},'_all\',Pair{nfolder,1},'Train',Pair{nfolder,2},'Test',CellScreen{nCellScreen},'\'];
%         load([result_path,'\',GroupName{nregion},'_all.mat'],'performanceRER')
%         Aved2(nregion,1,nfolder)=performanceRER(1,1);
%         Aved2(nregion,2,nfolder)=performanceRER(4,1);
%     end
% end
% Avem=mean(Aved2,3);


%% confusion matrix
OutPath=['E:\PFC and MTL\Results\Review\decoding\ByOcue\'];
GroupName = {'lPFC','mPFC','HPC'};
CellScreen = {'','_3','_correct','_correct_3'};
nCellScreen = 3;
%% all possible pairs
Pair={'-90','0';'-90','90';'0','-90';'0','90';'90','0';'90','-90'};
Decoder=[1 4];
DecoderName={'Item-location','Target-location'};
for nsub=2
     set(gcf,'OuterPosition',[992 248 1236 1075])
    for nregion = 1:3
        decodedLabels=[];trueLabels=[];
        for nfolder = 1:length(Pair(:,1))
            result_path = ['E:\PFC and MTL\Results\Review\decoding\ByOcue\Bar_OP\',GroupName{nregion},'_all\',Pair{nfolder,1},'Train',Pair{nfolder,2},'Test',CellScreen{nCellScreen},'\'];
            load([result_path,'\',GroupName{nregion},'_all.mat'],'PredictFit')
            decodedLabels=cat(1,decodedLabels,PredictFit{Decoder(nsub),1,1});
            trueLabels=cat(1,trueLabels,PredictFit{Decoder(nsub),1,2}(:,Decoder(nsub)));
        end
        subplot(1,3,nregion)
        tool_ConfusionMatrix_v2(decodedLabels,trueLabels)
        customMap=tool_BiColorMap(0.25);
        colormap(customMap)
        title([DecoderName{nsub},' ',GroupName{nregion}],'FontWeight','normal')
        set(gca,'FontSize',18)
    end
    pause
    FilePath=[OutPath,'/',DecoderName{nsub},'_BiColor.jpg'];
    saveas(gcf,FilePath)
    close all
end


%% confusion matrix for figure 4
OutPath=['E:\PFC and MTL\Results\Review\decoding\MixedDecoding\'];
GroupName = {'lPFC','mPFC','HPC'};
CellScreen = {'','_3','_correct','_correct_3'};
nCellScreen = 4;
%% all possible pairs
Pair={'-90','0';'-90','90';'0','-90';'0','90';'90','0';'90','-90'};
Decoder=[1 3 4];
DecoderName={'Item-location','Context-cue','Target-location'};
CutOff=[0.25 0.33 0.25];
for nsub=3
    figure()
    set(gcf,'OuterPosition',[992 248 1236 1075])
    for nregion = 1:3
            result_path = ['E:\PFC and MTL\Results\Review\decoding\MixedDecoding\OPMixed\',GroupName{nregion},'_all\',CellScreen{nCellScreen},'\'];
            load([result_path,'\',GroupName{nregion},'_all.mat'],'PredictFit')
            subplot(1,3,nregion)
            decodedLabels=PredictFit{Decoder(nsub),1,1};
            trueLabels=PredictFit{Decoder(nsub),1,2};
            tool_ConfusionMatrix_v2(decodedLabels,trueLabels)
            customMap=tool_BiColorMap(CutOff(nsub));
            colormap(customMap)
            colorbar
            title([DecoderName{nsub},' ',GroupName{nregion}],'FontWeight','normal')
            set(gca,'FontSize',18)
    end
    FilePath=[OutPath,'/',DecoderName{nsub},'_BiColor.jpg'];
    saveas(gcf,FilePath)
    close all
end

