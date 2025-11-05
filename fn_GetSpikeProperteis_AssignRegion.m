function fn_GetSpikeProperteis_AssignRegion(InPath,OutPath)

%% load spike information
load([InPath,'/SpikeProperties_v2_full.mat'],'CellInfo','CellName')

% concatenate all regions and animals
CellInfoAll=[];CellNameAll=[];
for nsec=1:length(CellName)
    CellInfoAll=cat(1,CellInfoAll,CellInfo{nsec});
    CellNameAll=cat(1,CellNameAll,CellName{nsec});
end

RegionName={'lPFC','mPFC','HPC'};
CellList=cell(3,1);
for nregion=1:3
    CellList{nregion}=importdata(['E:\PFC and MTL\cell name\',RegionName{nregion},'_all.txt']);
end

FinalCellName=cell(10,1);
FinalCellInfo=nan(10,4);
ncell=0;
for nsource=1:length(CellNameAll)
    for nregion=1:3
        if ismember(CellNameAll{nsource},CellList{nregion})
            ncell=ncell+1;
            FinalCellName{ncell}=CellNameAll{nsource};
            FinalCellInfo(ncell,1:3)=CellInfoAll(nsource,1:3);
            FinalCellInfo(ncell,4)=nregion;
        end
    end
end

% comparison
pair=[1 2;1 3;2 3];
for ncon=1:3
    for npair=1:3
        p(ncon,npair)=ranksum(FinalCellInfo(FinalCellInfo(:,4)==pair(npair,1),ncon),FinalCellInfo(FinalCellInfo(:,4)==pair(npair,2),ncon));
        Aved(ncon,npair,1)=nanmean(FinalCellInfo(FinalCellInfo(:,4)==pair(npair,1),ncon));
        Aved(ncon,npair,2)=nanmean(FinalCellInfo(FinalCellInfo(:,4)==pair(npair,2),ncon));
    end
end

save([OutPath,'/SpikeProperties_RegionAssigned.mat'],'FinalCellInfo','FinalCellName')
LogFinalCellInfo(:,1)=log10(FinalCellInfo(:,1));
LogFinalCellInfo(:,2)=log10(FinalCellInfo(:,2));
LogFinalCellInfo(:,3:4)=FinalCellInfo(:,3:4);


%% plot log

load([OutPath,'/SpikeProperties_RegionAssigned.mat'],'FinalCellInfo','FinalCellName')
LogFinalCellInfo(:,1)=log10(FinalCellInfo(:,1));
LogFinalCellInfo(:,2)=log10(FinalCellInfo(:,2));
LogFinalCellInfo(:,3:4)=FinalCellInfo(:,3:4);

perppair=[1 2;1 3;2 3];
ColorSet={[157,215,157]/255;[253,200,151]/255;[194,178,214]/255};
ColorSet_v2={[0, 0.6902, 0.3137];[0.9294, 0.4902, 0.1922];[0.4392, 0.1882, 0.6275]};
AxisName={'Firing rate (Hz)';'Burst fraction';'Spike width (ms)'};
XTICKS={[-2 -1 0 1 2];[-3 -2 -1 0 1];[0 0.2 0.4 0.6 0.8 1]};
XTICKLABEL={[0.01 0.1 1 10 100];[0.001,0.01 0.1 1 10];[0 0.2 0.4 0.6 0.8 1]};
% LIM={[-2 2];[-2 1];[0.2 1]};
ColorSet2={[0,1,0];[1 1 0];[1,0,1]};

for nsub=1:3
    subplot(1,3,nsub)
    for nregion=1:3
        indx=LogFinalCellInfo(:,4)==nregion;
        X=LogFinalCellInfo(indx,perppair(nsub,1));
        Y=LogFinalCellInfo(indx,perppair(nsub,2));
        scatter(X,Y,40,ColorSet{nregion},'filled','LineWidth',1,'MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6)
        hold on
        xlabel(AxisName{perppair(nsub,1)})
        ylabel(AxisName{perppair(nsub,2)})
        xticks(XTICKS{perppair(nsub,1)})
        xticklabels(XTICKLABEL{perppair(nsub,1)})
        yticks(XTICKS{perppair(nsub,2)})
        yticklabels(XTICKLABEL{perppair(nsub,2)})
        % xlim(LIM{perppair(nsub,1)})
        %  ylim(LIM{perppair(nsub,2)})
    end
    for nregion=1:3
        indx=LogFinalCellInfo(:,4)==nregion;
        X=median(LogFinalCellInfo(indx,perppair(nsub,1)),'omitnan');
        Y=median(LogFinalCellInfo(indx,perppair(nsub,2)),'omitnan');
        scatter(X,Y,80,ColorSet_v2{nregion},'hexagram','filled','LineWidth',0.75,'MarkerEdgeColor','k');
        hold on
    end
    legend lPFC mPFC HPC
    legend('Box','off');
    set(gca,'FontSize',18)
    set(gca, 'LineWidth', 2);
end
pause
FilePath=[OutPath,'/SpikeProperties_RegionAssigned_Log.jpg'];
saveas(gcf,FilePath)



%% plot raw
close all
perppair=[1 2;1 3;2 3];
ColorSet={[157,215,157]/255;[253,200,151]/255;[194,178,214]/255};
AxisName={'Firing rate (Hz)';'Burst fraction';'Spike width (ms)'};

for nsub=1:3
    subplot(1,3,nsub)
    for nregion=1:3
        indx=FinalCellInfo(:,4)==nregion;
        X=FinalCellInfo(indx,perppair(nsub,1));
        Y=FinalCellInfo(indx,perppair(nsub,2));
        scatter(X,Y,30,ColorSet{nregion},'filled')
        hold on
        xlabel(AxisName{perppair(nsub,1)})
        ylabel(AxisName{perppair(nsub,2)})
        legend lPFC mPFC HPC
        set(gca,'FontSize',18)
    end
end
FilePath=[OutPath,'/SpikeProperties_RegionAssigned.jpg'];
saveas(gcf,FilePath)