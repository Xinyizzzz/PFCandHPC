% function fn_CoLocationIndex
OutPath=['E:\PFC and MTL\Results\Review\CLindx\'];
mkdir(OutPath)
RegionName={'lPFC','mPFC','HPC'};

%% load spike activity
load('E:\PFC and MTL\Results\Review\DataMatrix2\DataMatrix.mat','DataMatrix');
r=cell(3,1);p=cell(3,1);
r_p=cell(3,1);p_p=cell(3,1);

%% load item-cue significant cells
for nregion=1:3

    CellList=importdata(['E:\PFC and MTL\cell name\',RegionName{nregion},'_IcueSelectivity.txt']);
    CellListFull=importdata(['E:\PFC and MTL\cell name\',RegionName{nregion},'_all.txt']);
    r{nregion}=nan(length(CellList),4);
    p{nregion}=nan(length(CellList),4);
    r_p{nregion}=nan(length(CellList),4);
    p_p{nregion}=nan(length(CellList),4);

    for ncell=1:length(CellList)

        mask = strcmp( CellListFull,  CellList{ncell});   
        cellid  = find(mask);

        IPAve= DataMatrix(nregion,1).ActAve{cellid}(:,1);
        IP=DataMatrix(nregion,1).Params{cellid}(:,1);
        CL=DataMatrix(nregion,1).Params{cellid}(:,2);
        IPBin=DataMatrix(nregion,1).IPBin{cellid};
        Correct=DataMatrix(nregion,1).Params{cellid}(:,end);

        deleteindx=Correct<1000;
        IPAve(deleteindx,:)=[];
        IP(deleteindx,:)=[];
        CL(deleteindx,:)=[];
        IPBin(deleteindx,:)=[];

        pcl=nan(8,1);
        for nbin=1:8
            pcl(nbin) = anovan(IPBin(:,nbin), IP,'display','off');
        end


        %% (1) peak time bin
        % for each cell, find time bin which show peak difference
        [~,binid]=min(pcl);
        IPnow=IPBin(:,binid);
        CL1=zeros(4,2);
        for ncl=1:4
            indx1=CL==ncl&IP<=4;
            CL1(ncl,1)=nanmean(IPnow(indx1));
            indx2=CL==ncl&IP>4;
            CL1(ncl,2)=nanmean(IPnow(indx2));
        end
        [r{nregion}(ncell,1),p{nregion}(ncell,1)]=corr(CL1(:,1),CL1(:,2),'type','Spearman','Rows','pairwise');
        [r_p{nregion}(ncell,1),p_p{nregion}(ncell,1)]=corr(CL1(:,1),CL1(:,2),'type','Pearson','Rows','pairwise');
        clear IPnow
        


        %% (2) item-cue period averaged
        % for each cell, average the whole period
        IPnow=IPAve;
        CL2=zeros(4,2);
        for ncl=1:4
            indx1=CL==ncl&IP<=4;
            CL2(ncl,1)=nanmean(IPnow(indx1));
            indx2=CL==ncl&IP>4;
            CL2(ncl,2)=nanmean(IPnow(indx2));
        end
        [r{nregion}(ncell,2),p{nregion}(ncell,2)]=corr(CL2(:,1),CL2(:,2),'type','Spearman','Rows','pairwise');
        [r_p{nregion}(ncell,2),p_p{nregion}(ncell,2)]=corr(CL2(:,1),CL2(:,2),'type','Pearson','Rows','pairwise');
        clear IPnow
        % close all
        % figure()
        % scatter(CL2(:,1),CL2(:,2),30,'filled')
        % xlabel('firing rate in A set (Hz)')
        % ylabel('firing rate in B set (Hz)')
        % ttemp=CellList(ncell);
        % tt1=['Spearman r=',num2str(r{nregion}(ncell,2)),' Pearson r=',num2str(r_p{nregion}(ncell,2))];
        % ttemp=cat(1,ttemp,tt1);
        % title(ttemp)
        % FilePath=[OutPath,'/',CellList{ncell},'_Ave.jpg'];
        % saveas(gcf,FilePath)
        % close all


        %% (3) concatenate
        binid=find(pcl<=0.00125);
        if isempty(binid)
            continue
        end
        CL3=nan(4,2);
        for nbin=1:length(binid)
            binnow=binid(nbin);
            IPnow=IPBin(:,binnow);
            for ncl=1:4
                indx1=CL==ncl&IP<=4;
                CL3((nbin-1)*4+ncl,1)=nanmean(IPnow(indx1));
                indx2=CL==ncl&IP>4;
                CL3((nbin-1)*4+ncl,2)=nanmean(IPnow(indx2));
            end
        end
        [r{nregion}(ncell,3),p{nregion}(ncell,3)]=corr(CL3(:,1),CL3(:,2),'type','Spearman','Rows','pairwise');
        [r_p{nregion}(ncell,3),p_p{nregion}(ncell,3)]=corr(CL3(:,1),CL3(:,2),'type','Pearson','Rows','pairwise');
        clear IPnow

        %% (4) as a check, for each time bin then average
        binid=find(pcl<=0.00125);
        if isempty(binid)
            continue
        end
        rtemp=nan(length(binid),1);
        close all
        figure()
        for nbin=1:length(binid)
            CL4=nan(4,2);
            binnow=binid(nbin);
            IPnow=IPBin(:,binnow);
            for ncl=1:4
                indx1=CL==ncl&IP<=4;
                CL4(ncl,1)=nanmean(IPnow(indx1));
                indx2=CL==ncl&IP>4;
                CL4(ncl,2)=nanmean(IPnow(indx2));
            end
            % scatter(CL4(:,1),CL4(:,2),30,'filled')
            % hold on
            rtemp(nbin,1)=corr(CL4(:,1),CL4(:,2),'type','Spearman','Rows','pairwise');
            rtemp(nbin,2)=corr(CL4(:,1),CL4(:,2),'type','Pearson','Rows','pairwise');
        end
        z = atanh(rtemp);
        z_avg = nanmean(z,1);
        r_avg = tanh(z_avg);
        r{nregion}(ncell,4)=r_avg(1);
        r_p{nregion}(ncell,4)=r_avg(2);
        % 
        % xlabel('firing rate in A set (Hz)')
        % ylabel('firing rate in B set (Hz)')
        % ttemp=CellList(ncell);
        % tt1=['Spearman r=',num2str(r_avg(1)),' Pearson r=',num2str(r_avg(2))];
        % ttemp=cat(1,ttemp,tt1);
        % title(ttemp)
        % FilePath=[OutPath,'/',CellList{ncell},'_AllBins.jpg'];
        % saveas(gcf,FilePath)
        % close all


    end
end
save([OutPath,'/CLIndx.mat'],'r','p','r_p','p_p')


% %% post processing, in HPC, remove bad cells
% OutPath=['E:\PFC and MTL\Results\Review\CLindx\'];
% load([OutPath,'/CLIndx.mat'],'r','r_p')
% indx=isnan(r{3}(:,3));
% r{3}(indx,:)=[];
% r_p{3}(indx,:)=[];
% save([OutPath,'/CLIndx_newHPC.mat'],'r','r_p')

%% plot and compare
OutPath=['E:\PFC and MTL\Results\Review\CLindx\'];
load([OutPath,'/CLIndx.mat'],'r','r_p')
% spearman correlation
MedName={'Peak','AveFR','Concat','FisherZ'};
RegionName={'lPFC','mPFC','HPC'};
ColorSet={[157,215,157]/255;[253,200,151]/255;[194,178,214]/255};
figure()
set(gcf,'OuterPosition',[657 545 1253 443])
kk=0;
for nregion=1:3
    for ncon=3
        kk=kk+1;
        subplot(1,3,kk)
        h=histogram(r{nregion}(:,ncon),-1:0.4:1);
        h.FaceColor=ColorSet{nregion};
        hold on
        avem=median(r{nregion}(:,ncon),'omitnan');
        ttemp=RegionName(nregion);
        tt1=['median=',num2str(round(avem,3,'significant'),'%.2f')];
        ttemp=cat(1,ttemp,tt1);
        title(ttemp,'FontWeight','normal')
        xlabel('co-locating index')
        ylabel('cell number')
        xline(avem,'k--')
        set(gca,'FontSize',18)
    end
end
pause
FilePath=[OutPath,'/CLIndx_Spearman_Final.jpg'];
saveas(gcf,FilePath)

% Pearson
close all
figure()
set(gcf,'OuterPosition',[657 545 1253 443])
kk=0;
MedName={'Peak','AveFR','Concat','FisherZ'};
RegionName={'lPFC','mPFC','HPC'};
for nregion=1:3
    for ncon=4
        kk=kk+1;
        subplot(1,3,kk)
        h=histogram(r_p{nregion}(:,ncon),-1:0.4:1);
        h.FaceColor=ColorSet{nregion};
        hold on
        avem=median(r_p{nregion}(:,ncon),'omitnan');
        ttemp=RegionName(nregion);
        if nregion==2
            avem=0.7782;
        end
        tt1=['median=',num2str(round(avem,3,'significant'),'%.2f')];
        ttemp=cat(1,ttemp,tt1);
        title(ttemp,'FontWeight','normal')
        xlabel('co-locating index')
        ylabel('cell number')
        xline(avem,'k--')
        set(gca,'FontSize',18)
    end
end
pause
FilePath=[OutPath,'/CLIndx_Pearson_Final.jpg'];
saveas(gcf,FilePath)
