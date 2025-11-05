function fn_EyePos_Correlation_Recenter
InPath=['E:\PFC and MTL\Results\Review\EyeDataMatrix_noShift\'];
OutPath=['E:\PFC and MTL\Results\Review\EyePosition_noShift\Summary\Raw\Corrected\Recentered\'];
mkdir(OutPath)

load([InPath,'/DataMatrix_wEye.mat'],'DataMatrix');
Summary=cell(3,3,2,4);


for nregion=1:size(DataMatrix,1)
    if nregion<=2
        ParID{1}=[1 4 3 2 5 8 7 6];
    else
        ParID{1}=[1 2 3 4 5 6 7 8];
    end
    if nregion<=2
        ParID{2}=[1 4 3 2];
    else
        ParID{2}=[1 2 3 4];
    end
    ParID{3}=[-90 0 90];
    if nregion<=2
        ParID{4}=[1 4 3 2];
    else
        ParID{4}=[1 3 5 7];
    end

    for ndim=1:2
        for nperiod=1:3
            for npar=1:4
                Summary{nregion,nperiod,ndim,npar}=nan(size(DataMatrix(nregion).Params,2),size(ParID{npar},2));
            end
        end
    end

    for ndim=1:2
        for nperiod=1:3
            % Params=[Icue,CL,Ocue,Tar,Choice,Correct];
            % DataMatrix(nregion,1).EyeAve{ncell,1}=EyeAve_H;
            % DataMatrix(nregion,1).EyeAve{ncell,2}=EyeAve_V;
            % EyeAve_H=squeeze(eye_h_fin(:,32)) IP;33 OP

            for ncell=1:size(DataMatrix(nregion).Params,2)

                Params=DataMatrix(nregion).Params{ncell};

                if  length(DataMatrix(nregion,1).EyeAve{ncell,ndim})<20
                    continue
                end

                Correct=Params(:,end);
                Ocue=Params(:,3);
                indx=Correct>999&abs(Ocue)~=45;
                Ave=DataMatrix(nregion,1).EyeAve{ncell,ndim}(indx,nperiod);
                % if nperiod==1
                %     Zve=zscore(Ave);
                % else
                %     Zve=Ave-DataMatrix(nregion,1).EyeAve{ncell,ndim}(indx,nperiod-1);
                % end
                Zve=Ave-nanmean(Ave);
                % Zve=zscore(Ave);
                % Zve=Ave;
                Params_filtered=Params(indx,:);

                for npar=1:4
                    for nid=1:size(ParID{npar},2)
                        % Summary{nregion,nperiod,ndim,npar}(ncell,nid)=nanmean(Zve(Params_filtered(:,npar)==ParID{npar}(nid)),1);
                        if nregion>2
                            Summary{nregion,nperiod,ndim,npar}(ncell,nid)=nanmean(Zve(Params_filtered(:,npar)==ParID{npar}(nid)),1);
                        else
                            Summary{nregion,nperiod,ndim,npar}(ncell,nid)=-1*nanmean(Zve(Params_filtered(:,npar)==ParID{npar}(nid)),1);
                        end
                    end
                end
            end
        end
    end
end

save([OutPath,'/PerformanceSummary.mat'],'Summary')



%% resort by animal
OutPath=['E:\PFC and MTL\Results\Review\EyePosition_noShift\Summary\Raw\Corrected\Recentered\'];
load([OutPath,'/PerformanceSummary.mat'],'Summary')
RegionName={'lPFC','mPFC','HPC'};
AnimalName={'Care','Galeno','Billy'};
SummaryEachAnimal=cell(3,3,2,4);
for nregion=1:3
    CellList=importdata(['E:\PFC and MTL\cell name\',RegionName{nregion},'_eye.txt']);
    for ncell=1:length(CellList)
        for nanimal=1:3
            if strcmp(CellList{ncell}(1),AnimalName{nanimal}(1))
                Animal=nanimal;
                break
            end
        end
        %Summary{nregion,nperiod,ndim,npar}
        for nperiod=1:3
            for ndim=1:2
                for npar=1:4
                    SummaryEachAnimal{Animal,nperiod,ndim,npar}=cat(1,SummaryEachAnimal{Animal,nperiod,ndim,npar},Summary{nregion,nperiod,ndim,npar}(ncell,:));
                end
            end
        end
    end
end
save([OutPath,'/PerformanceSummary_EachAnimal.mat'],'SummaryEachAnimal')


%% plot each animal to confirm location
clear all
clc
OutPath=['E:\PFC and MTL\Results\Review\EyePosition_noShift\Summary\Raw\Corrected\Recentered\'];
load([OutPath,'/PerformanceSummary_EachAnimal.mat'],'SummaryEachAnimal')
Summary=SummaryEachAnimal;
RegionName={'Monkey C','Monkey G','Monkey B'};
PeriodName={'item-cue','context-cue','choice'};
ParName={'Item identity','Item-location','Context-cue','Target-location'};
DimName={'H','V'};
XTICKLABEL{1}={'I-A','II-A','III-A','IV-A','I-B','II-B','III-B','IV-B'};
XTICKLABEL{2}={'I','II','III','IV'};
XTICKLABEL{3}={'-90','0','90'};
XTICKLABEL{4}={'TR','BR','BL','TL'};
Sym={'o','^','square'};
RGB=orderedcolors('gem');
RGB(5:8,:)=RGB(1:4,:);
p=nan(2,4,2,4);
for nperiod=3
    for npar=4
        MonkeyName=cell(1,1);
        Condition=cell(1,1);


        for nregion=1:3
            %Summary{nregion,nperiod,ndim,npar}(ncell,nid)
            figure()
            set(gcf,'OuterPosition',[992 843 430 480])
            AA=cell(size(XTICKLABEL{npar},2),1);
            for nid=1:size(XTICKLABEL{npar},2)
                AA{nid}(:,1)=Summary{nregion,nperiod,1,npar}(:,nid);
                AA{nid}(:,2)=Summary{nregion,nperiod,2,npar}(:,nid);

                Avex=nanmean(AA{nid}(:,1));
                Avey=nanmean(AA{nid}(:,2));
                stdx=std(AA{nid}(:,1),'omitnan');
                stdy=std(AA{nid}(:,2),'omitnan');
                hold on

                errorbar(Avex,Avey,stdx,'horizontal','Color',RGB(nid,:))
                errorbar(Avex,Avey,stdy,'Color',RGB(nid,:))
                scatter(Avex,Avey,60,RGB(nid,:),'filled','Marker',Sym{nregion})

                MonkeyName=cat(1,MonkeyName,RegionName{nregion});
                Condition=cat(1,Condition,XTICKLABEL{npar}{nid});
            end

            xlabel(['Horizontal position (dva)'])
            ylabel(['Vertical position (dva)'])
            ttemp1={[PeriodName{nperiod},' period ']};
            tt1=ParName{npar};
            ttemp=cat(1,tt1,ttemp1);

            xline(0,'k--')
            yline(0,'k--')
            % xlim([-1 1])
            % ylim([-1 1])
            % xticks([-1 -0.5 0 0.5 1])
            % yticks([-1 -0.5 0 0.5 1])
             xlim([-10 10])
            ylim([-10 10])
            xticks([-10 -5 0 5 10])
            yticks([-10 -5 0 5 10])

            set(gca,'FontSize',18)
            title(ttemp,'FontSize',17,'FontWeight','normal')


            pause
            ttsave=[RegionName{nregion},' ',PeriodName{nperiod},'-period ',ParName{npar}];
            FilePath=[OutPath,'/Plot/ZoomInIn_scatter_',ttsave,'.jpg'];
            saveas(gcf,FilePath)
            close all
        end
    end
end


%% cosine similarity between co-locating items
clear all
clc
close all
OutPath=['E:\PFC and MTL\Results\Review\EyePosition_noShift\Summary\Raw\Corrected\Recentered\'];
load([OutPath,'/PerformanceSummary_EachAnimal.mat'],'SummaryEachAnimal')
Summary=SummaryEachAnimal;
%Summary{nregion,nperiod,ndim,npar}(ncell,nid)
CosSim=cell(3,3);
InnerP=cell(3,3);
Ave=nan(3,3,4);
for nregion=1:3
    for nperiod=1:3
        CosSim{nregion,nperiod}=nan(size(Summary{nregion,1,1,1},1),4);
        for npar=1
            for ncell=1:size(Summary{nregion,nperiod,1,npar},1)
                for nid=1:4
                    clear A B
                    A(1,1)=Summary{nregion,nperiod,1,npar}(ncell,nid);
                    A(1,2)=Summary{nregion,nperiod,2,npar}(ncell,nid);
                    B(1,1)=Summary{nregion,nperiod,1,npar}(ncell,nid+4);
                    B(1,2)=Summary{nregion,nperiod,2,npar}(ncell,nid+4);
                    CosSim{nregion,nperiod}(ncell,nid)=dot(A,B)/(norm(A)*norm(B));
                    InnerP{nregion,nperiod}(ncell,nid)=dot(A,B);
                end
            end
        end
        Ave(nregion,nperiod,1:4)=nanmean(CosSim{nregion,nperiod});
    end
end
save([OutPath,'/CosSimilarity.mat'],'CosSim','Ave','InnerP')




%% plot and statistics
% within each animal
DimName={'Horizontal','Vertical'};
PeriodName={'item-cue period','context-cue period','choice fixation'};
for nperiod=1:3
    p=nan(4,1);
    for nregion=1:3
        Ave=nanmean(sum(CosSim{nregion,nperiod},2))/4;
        STD=std(sum(CosSim{nregion,nperiod},2),'omitnan')/4;
        bar(nregion,Ave)
        hold on
        errorbar(nregion,Ave,STD)
        p(nregion,1)=signrank(sum(CosSim{nregion,nperiod},2),0,'tail','right');
    end
    xticks([1 2 3])
    xticklabels({'Care','Galeno','Billy'})
    % ylim([-1 1])
    ttemp=[PeriodName{nperiod}];
    tt1=['pcare=',num2str(p(1)),' pg=',num2str(p(2)),' pb=',num2str(p(3))];
    tt1=cat(1,{ttemp},tt1);
    title(tt1);
    pause
    FilePath=[OutPath,'/ZCosSimilaritySum_',ttemp,'.jpg'];
    saveas(gcf,FilePath)
    close all
end

%% correlation separately for horizontal and vertical location for co-locating items
OutPath=['E:\PFC and MTL\Results\Review\EyePosition_noShift\Summary\Raw\Corrected\Recentered\'];
load([OutPath,'/PerformanceSummary_EachAnimal.mat'],'SummaryEachAnimal')
Summary=SummaryEachAnimal;
%Summary{nregion,nperiod,ndim,npar}(ncell,nid)
P=cell(3,3);
Ave=nan(3,3,2);
for nregion=1:3
    for nperiod=1:3
        P{nregion,nperiod}=nan(size(Summary{nregion,1,1,1},1),2);
        for ndim=1:2
            for npar=1
                for ncell=1:size(Summary{nregion,nperiod,ndim,npar},1)
                    clear A B
                    A=Summary{nregion,nperiod,ndim,npar}(ncell,[1 2 3 4]);
                    B=Summary{nregion,nperiod,ndim,npar}(ncell,[5 6 7 8]);
                    A=A(:);
                    B=B(:);
                    P{nregion,nperiod}(ncell,ndim)=corr(A,B,'rows','pairwise');
                end
            end
        end
        Ave(nregion,nperiod,1:2)=nanmean(P{nregion,nperiod});
    end
end
save([OutPath,'/Correlation.mat'],'P','Ave')


%% plot and statistics
% within each animal
DimName={'Horizontal','Vertical'};
PeriodName={'item-cue period','context-cue period','choice fixation'};
for ndim=1:2
    for nperiod=1
        p=nan(4,1);
        for nregion=1:3
            Ave=nanmean(P{nregion,nperiod}(:,ndim));
            STD=std(P{nregion,nperiod}(:,ndim),'omitnan');
            bar(nregion,Ave)
            hold on
            errorbar(nregion,Ave,STD)
            p(nregion,1)=signrank(P{nregion,nperiod}(:,ndim),0);
        end
        xticks([1 2 3])
        xticklabels({'Care','Galeno','Billy'})
        % ylim([-1 1])
        ttemp=[DimName{ndim},' ',PeriodName{nperiod}];
        tt1=['pcare=',num2str(p(1)),' pg=',num2str(p(2)),' pb=',num2str(p(3))];
        tt1=cat(1,{ttemp},tt1);
        title(tt1);
        pause
        FilePath=[OutPath,'/ColocationItemCorr_',ttemp,'.jpg'];
        saveas(gcf,FilePath)
        close all

    end
end


%% cosine similarity between period

clear all
clc
close all
OutPath=['E:\PFC and MTL\Results\Review\EyePosition_noShift\Summary\Raw\Corrected\Recentered\'];
load([OutPath,'/PerformanceSummary_EachAnimal.mat'],'SummaryEachAnimal')
Summary=SummaryEachAnimal;
%Summary{nregion,nperiod,ndim,npar}(ncell,nid)
CosSim=cell(3,3);
InnerP=cell(3,3);
Ave=nan(3,3,5);
for nregion=1:3
    kk=0;
    for nperiod1=1:2
        for nperiod2=nperiod1+1:3
            kk=kk+1;
            CosSim{nregion,kk}=nan(size(Summary{nregion,1,1,1},1),4);
            for ncell=1:size(Summary{nregion,1,1,1},1)
                for nid=1:4
                    clear A B
                    if nperiod1==1
                        A(1,1)=Summary{nregion,nperiod1,1,2}(ncell,nid);
                        A(1,2)=Summary{nregion,nperiod1,2,2}(ncell,nid);
                    else
                        A(1,1)=Summary{nregion,nperiod1,1,4}(ncell,nid);
                        A(1,2)=Summary{nregion,nperiod1,2,4}(ncell,nid);
                    end
                    B(1,1)=Summary{nregion,nperiod2,1,4}(ncell,nid);
                    B(1,2)=Summary{nregion,nperiod2,2,4}(ncell,nid);
                    CosSim{nregion,kk}(ncell,nid)=dot(A,B)/(norm(A)*norm(B));
                    InnerP{nregion,kk}(ncell,nid)=dot(A,B);
                end
            end
            Ave(nregion,kk,1:4)=nanmean(CosSim{nregion,kk});
        end
    end
    kk=kk+1;
    CosSim{nregion,kk}=nan(size(Summary{nregion,1,1,1},1),4);
    for ncell=1:size(Summary{nregion,1,1,1},1)
        for nid=1:4
            clear A B
            A(1,1)=Summary{nregion,2,1,2}(ncell,nid);
            A(1,2)=Summary{nregion,2,2,2}(ncell,nid);
            B(1,1)=Summary{nregion,3,1,4}(ncell,nid);
            B(1,2)=Summary{nregion,3,2,4}(ncell,nid);
            CosSim{nregion,kk}(ncell,nid)=dot(A,B)/(norm(A)*norm(B));
            InnerP{nregion,kk}(ncell,nid)=dot(A,B);
        end
    end
    Ave(nregion,kk,1:4)=nanmean(CosSim{nregion,kk});
    kk=kk+1;
    CosSim{nregion,kk}=nan(size(Summary{nregion,1,1,1},1),4);
    for ncell=1:size(Summary{nregion,1,1,1},1)
        for nid=1:4
            clear A B
            A(1,1)=Summary{nregion,2,1,2}(ncell,nid);
            A(1,2)=Summary{nregion,2,2,2}(ncell,nid);
            B(1,1)=Summary{nregion,2,1,4}(ncell,nid);
            B(1,2)=Summary{nregion,2,2,4}(ncell,nid);
            CosSim{nregion,kk}(ncell,nid)=dot(A,B)/(norm(A)*norm(B));
            InnerP{nregion,kk}(ncell,nid)=dot(A,B);
        end
    end
    Ave(nregion,kk,1:4)=nanmean(CosSim{nregion,kk});
end
save([OutPath,'/3Period_CosSimilarity.mat'],'CosSim','Ave','InnerP')


%% plot and statistics
% within each animal
DimName={'Horizontal','Vertical'};
PeriodName={'Icue VS Ocue','Icue VS choice','Ocue VS choice', 'OcueCL VS choice','OcueCL VS OCueTar'};
for nperiod=4
    p=nan(4,1);
    for nregion=1:3
        P=[];
        % P=cat(1,CosSim{nregion,nperiod}(:,1),CosSim{nregion,nperiod}(:,2),CosSim{nregion,nperiod}(:,3),CosSim{nregion,nperiod}(:,4));
        P=nanmean(CosSim{nregion,nperiod},2);
        Ave=nanmean(P);
        STD=std(P,'omitnan');
        bar(nregion,Ave)
        hold on
        errorbar(nregion,Ave,STD)
        p(nregion,1)=signrank(P,0,'tail','right');
    end
    xticks([1 2 3])
    xticklabels({'Care','Galeno','Billy'})
    % ylim([-1 1])
    ttemp=[PeriodName{nperiod}];
    tt1=['pcare=',num2str(p(1)),' pg=',num2str(p(2)),' pb=',num2str(p(3))];
    tt1=cat(1,{ttemp},tt1);
    title(tt1);
    pause
    FilePath=[OutPath,'/DCosSimilaritySum_',ttemp,'.jpg'];
    saveas(gcf,FilePath)
    close all
end


%% correlation between different period, separately for horizontal location and vertical location
OutPath=['E:\PFC and MTL\Results\Review\EyePosition_noShift\Summary\Raw\Corrected\Recentered\'];
load([OutPath,'/PerformanceSummary_EachAnimal.mat'],'SummaryEachAnimal')
Summary=SummaryEachAnimal;
%Summary{nregion,nperiod,ndim,npar}(ncell,nid)
P=cell(3,3);
Ave=nan(3,3,2);
for nregion=1:3
    kk=0;
    for nperiod1=1:2
        for nperiod2=nperiod1+1:3
            kk=kk+1;
            P{nregion,kk}=nan(size(Summary{nregion,1,1,1},1),2);
            for ndim=1:2
                for ncell=1:size(Summary{nregion,1,ndim,1},1)
                    clear A B
                    if nperiod1==1
                        A=Summary{nregion,nperiod1,ndim,2}(ncell,[1 2 3 4]);
                    else
                        A=Summary{nregion,nperiod1,ndim,4}(ncell,[1 2 3 4]);
                    end
                    B=Summary{nregion,nperiod2,ndim,4}(ncell,[1 2 3 4]);
                    A=A(:);
                    B=B(:);
                    P{nregion,kk}(ncell,ndim)=corr(A,B,'rows','pairwise');
                end
            end
        end
        Ave(nregion,kk,1:2)=nanmean(P{nregion,kk});
    end
end
save([OutPath,'/Correlation_3Period.mat'],'P','Ave')

%% plot and statistics
% within each animal
DimName={'Horizontal','Vertical'};
PeriodName={'Icue VS Ocue','Icue VS choice','Ocue VS choice'};
for ndim=1:2
    for nperiod=3
        p=nan(4,1);
        for nregion=1:3
            Ave=nanmean(P{nregion,nperiod}(:,ndim));
            STD=std(P{nregion,nperiod}(:,ndim),'omitnan');
            bar(nregion,Ave)
            hold on
            errorbar(nregion,Ave,STD)
            p(nregion,1)=signrank(P{nregion,nperiod}(:,ndim),0);
        end
        xticks([1 2 3])
        xticklabels({'Care','Galeno','Billy'})
        % ylim([-1 1])
        ttemp=[DimName{ndim},' ',PeriodName{nperiod}];
        tt1=['pcare=',num2str(p(1)),' pg=',num2str(p(2)),' pb=',num2str(p(3))];
        tt1=cat(1,{ttemp},tt1);
        title(tt1);
        pause
        FilePath=[OutPath,'/3Period_',ttemp,'.jpg'];
        saveas(gcf,FilePath)
        close all

    end
end


%% for each parameter, does the cosine similarity between each category (e.g., between co-location I and co-location II) as high as the choice-fixation period?
clear all
clc
close all
OutPath=['E:\PFC and MTL\Results\Review\EyePosition_noShift\Summary\Raw\Corrected\Recentered\'];
load([OutPath,'/PerformanceSummary_EachAnimal.mat'],'SummaryEachAnimal')
Summary=SummaryEachAnimal;
%Summary{nregion,nperiod,ndim,npar}(ncell,nid)
CosSim=cell(3,4,3);
for nregion=1:3
    for nperiod=1:3
        for npar=1:4
            CosSim{nregion,npar,nperiod}=nan(size(Summary{nregion,nperiod,1,npar},1),(size(Summary{nregion,nperiod,1,npar},2)-1)*size(Summary{nregion,nperiod,1,npar},2)/2);
            kk=0;
            for nid1=1:size(Summary{nregion,1,1,1,1},2)-1
                for nid2=(nid1+1):size(Summary{nregion,nperiod,1,npar},2)
                    kk=kk+1;
                    for ncell=1:size(Summary{nregion,1,1,1},1)
                        clear A B
                        A(1,1)=Summary{nregion,nperiod,1,npar}(ncell,nid1);
                        A(1,2)=Summary{nregion,nperiod,2,npar}(ncell,nid1);
                        B(1,1)=Summary{nregion,nperiod,1,npar}(ncell,nid2);
                        B(1,2)=Summary{nregion,nperiod,2,npar}(ncell,nid2);
                        CosSim{nregion,npar,nperiod}(ncell,kk)=dot(A,B)/(norm(A)*norm(B));
                    end
                end
            end
        end
    end
end
save([OutPath,'/WithinPeriodCrossItem_CosSimilarity.mat'],'CosSim')


%%%%%%%%%%%%%%%%%%%% for each region %%%%%%%%%%%%%%%%%
%% cosine similarity between period

clear all
clc
close all
OutPath=['E:\PFC and MTL\Results\Review\EyePosition_noShift\Summary\Raw\Corrected\Recentered\'];
load([OutPath,'/PerformanceSummary.mat'],'Summary')
%Summary{nregion,nperiod,ndim,npar}(ncell,nid)
CosSim=cell(3,3);
InnerP=cell(3,3);
Ave=nan(3,3,4);
for nregion=1:3
    kk=0;
    for nperiod1=1:2
        for nperiod2=nperiod1+1:3
            kk=kk+1;
            CosSim{nregion,kk}=nan(size(Summary{nregion,1,1,1},1),4);
            for ncell=1:size(Summary{nregion,1,1,1},1)
                for nid=1:4
                    clear A B
                    if nperiod1==1
                        A(1,1)=Summary{nregion,nperiod1,1,2}(ncell,nid);
                        A(1,2)=Summary{nregion,nperiod1,2,2}(ncell,nid);
                    else
                        A(1,1)=Summary{nregion,nperiod1,1,4}(ncell,nid);
                        A(1,2)=Summary{nregion,nperiod1,2,4}(ncell,nid);
                    end
                    B(1,1)=Summary{nregion,nperiod2,1,4}(ncell,nid);
                    B(1,2)=Summary{nregion,nperiod2,2,4}(ncell,nid);
                    CosSim{nregion,kk}(ncell,nid)=dot(A,B)/(norm(A)*norm(B));
                    InnerP{nregion,kk}(ncell,nid)=dot(A,B);
                end
            end
            Ave(nregion,kk,1:4)=nanmean(CosSim{nregion,kk});
        end
    end
    kk=kk+1;
    CosSim{nregion,kk}=nan(size(Summary{nregion,1,1,1},1),4);
    for ncell=1:size(Summary{nregion,1,1,1},1)
        for nid=1:4
            clear A B
            A(1,1)=Summary{nregion,2,1,2}(ncell,nid);
            A(1,2)=Summary{nregion,2,2,2}(ncell,nid);
            B(1,1)=Summary{nregion,3,1,4}(ncell,nid);
            B(1,2)=Summary{nregion,3,2,4}(ncell,nid);
            CosSim{nregion,kk}(ncell,nid)=dot(A,B)/(norm(A)*norm(B));
            InnerP{nregion,kk}(ncell,nid)=dot(A,B);
        end
    end

    Ave(nregion,kk,1:4)=nanmean(CosSim{nregion,kk});
end
save([OutPath,'/3Period_CosSimilarity_Region.mat'],'CosSim','Ave','InnerP')


%% plot and statistics
% within each animal
DimName={'Horizontal','Vertical'};
PeriodName={'Icue VS Ocue','Icue VS choice','Ocue VS choice', 'OcueCL VS choice'};
for nperiod=1:4
    p=nan(4,1);
    for nregion=1:3
        P=[];
        % P=cat(1,CosSim{nregion,nperiod}(:,1),CosSim{nregion,nperiod}(:,2),CosSim{nregion,nperiod}(:,3),CosSim{nregion,nperiod}(:,4));
        P=nanmean(CosSim{nregion,nperiod},2);
        Ave=nanmean(P);
        STD=std(P,'omitnan');
        bar(nregion,Ave)
        hold on
        errorbar(nregion,Ave,STD)
        p(nregion,1)=signrank(P,0,'tail','right');
    end
    xticks([1 2 3])
    xticklabels({'lPFC','mPFC','HPC'})
    % ylim([-1 1])
    ttemp=[PeriodName{nperiod}];
    tt1=['pcare=',num2str(p(1)),' pg=',num2str(p(2)),' pb=',num2str(p(3))];
    tt1=cat(1,{ttemp},tt1);
    title(tt1);
    pause
    FilePath=[OutPath,'/Region_CosSimilaritySum_',ttemp,'.jpg'];
    saveas(gcf,FilePath)
    close all
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% shuffle %%%%%%%%%%%%%%%%%%%%%%
%% shuffle

InPath=['E:\PFC and MTL\Results\Review\EyeDataMatrix_noShift\'];
OutPath=['E:\PFC and MTL\Results\Review\EyePosition_noShift\Summary\Raw\Corrected\Recentered\'];
mkdir(OutPath)

load([InPath,'/DataMatrix_wEye.mat'],'DataMatrix');
PermutationAve=nan(3,4,4,100);

for nshuffle=1:100

    disp(nshuffle)

    Summary=cell(3,3,2,4);


    for nregion=1:size(DataMatrix,1)
        if nregion<=2
            ParID{1}=[1 4 3 2 5 8 7 6];
        else
            ParID{1}=[1 2 3 4 5 6 7 8];
        end
        if nregion<=2
            ParID{2}=[1 4 3 2];
        else
            ParID{2}=[1 2 3 4];
        end
        ParID{3}=[-90 0 90];
        if nregion<=2
            ParID{4}=[1 4 3 2];
        else
            ParID{4}=[1 3 5 7];
        end

        for ndim=1:2
            for nperiod=1:3
                for npar=1:4
                    Summary{nregion,nperiod,ndim,npar}=nan(size(DataMatrix(nregion).Params,2),size(ParID{npar},2));
                end
            end
        end

        for ndim=1:2
            for nperiod=1:3
                % Params=[Icue,CL,Ocue,Tar,Choice,Correct];
                % DataMatrix(nregion,1).EyeAve{ncell,1}=EyeAve_H;
                % DataMatrix(nregion,1).EyeAve{ncell,2}=EyeAve_V;
                % EyeAve_H=squeeze(eye_h_fin(:,32)) IP;33 OP

                for ncell=1:size(DataMatrix(nregion).Params,2)

                    Params=DataMatrix(nregion).Params{ncell};

                    if  length(DataMatrix(nregion,1).EyeAve{ncell,ndim})<20
                        continue
                    end

                    Correct=Params(:,end);
                    Ocue=Params(:,3);
                    indx=Correct>999&abs(Ocue)~=45;
                    Ave=DataMatrix(nregion,1).EyeAve{ncell,ndim}(indx,nperiod);
                    if nperiod==1
                        Zve=zscore(Ave);
                    else
                        Zve=Ave-DataMatrix(nregion,1).EyeAve{ncell,ndim}(indx,nperiod-1);
                    end
                    % Zve=Ave-nanmean(Ave);
                    % Zve=zscore(Ave);
                    Params_filtered=Params(indx,:);
                    Zve=Zve(randperm(length(Zve)));

                    for npar=1:4
                        for nid=1:size(ParID{npar},2)
                            % Summary{nregion,nperiod,ndim,npar}(ncell,nid)=nanmean(Zve(Params_filtered(:,npar)==ParID{npar}(nid)),1);
                            if nregion>2
                                Summary{nregion,nperiod,ndim,npar}(ncell,nid)=nanmean(Zve(Params_filtered(:,npar)==ParID{npar}(nid)),1);
                            else
                                Summary{nregion,nperiod,ndim,npar}(ncell,nid)=-1*nanmean(Zve(Params_filtered(:,npar)==ParID{npar}(nid)),1);
                            end
                        end
                    end
                end
            end
        end
    end

    % resort by animal
    RegionName={'lPFC','mPFC','HPC'};
    AnimalName={'Care','Galeno','Billy'};
    SummaryEachAnimal=cell(3,3,2,4);
    for nregion=1:3
        CellList=importdata(['E:\PFC and MTL\cell name\',RegionName{nregion},'_eye.txt']);
        for ncell=1:length(CellList)
            for nanimal=1:3
                if strcmp(CellList{ncell}(1),AnimalName{nanimal}(1))
                    Animal=nanimal;
                    break
                end
            end
            %Summary{nregion,nperiod,ndim,npar}
            for nperiod=1:3
                for ndim=1:2
                    for npar=1:4
                        SummaryEachAnimal{Animal,nperiod,ndim,npar}=cat(1,SummaryEachAnimal{Animal,nperiod,ndim,npar},Summary{nregion,nperiod,ndim,npar}(ncell,:));
                    end
                end
            end
        end
    end

    % cosine similarity
    Summary=SummaryEachAnimal;
    %Summary{nregion,nperiod,ndim,npar}(ncell,nid)
    CosSim=cell(3,3);
    Ave=nan(3,3,4);
    for nregion=1:3
        kk=0;
        for nperiod1=1:2
            for nperiod2=nperiod1+1:3
                kk=kk+1;
                CosSim{nregion,kk}=nan(size(Summary{nregion,1,1,1},1),4);
                for ncell=1:size(Summary{nregion,1,1,1},1)
                    for nid=1:4
                        clear A B
                        if nperiod1==1
                            A(1,1)=Summary{nregion,nperiod1,1,2}(ncell,nid);
                            A(1,2)=Summary{nregion,nperiod1,2,2}(ncell,nid);
                        else
                            A(1,1)=Summary{nregion,nperiod1,1,4}(ncell,nid);
                            A(1,2)=Summary{nregion,nperiod1,2,4}(ncell,nid);
                        end
                        B(1,1)=Summary{nregion,nperiod2,1,4}(ncell,nid);
                        B(1,2)=Summary{nregion,nperiod2,2,4}(ncell,nid);
                        CosSim{nregion,kk}(ncell,nid)=dot(A,B)/(norm(A)*norm(B));
                    end
                end
                Ave(nregion,kk,1:4)=nanmean(CosSim{nregion,kk});
            end
        end
        kk=kk+1;
        CosSim{nregion,kk}=nan(size(Summary{nregion,1,1,1},1),4);
        for ncell=1:size(Summary{nregion,1,1,1},1)
            for nid=1:4
                clear A B
                A(1,1)=Summary{nregion,2,1,2}(ncell,nid);
                A(1,2)=Summary{nregion,2,2,2}(ncell,nid);
                B(1,1)=Summary{nregion,3,1,4}(ncell,nid);
                B(1,2)=Summary{nregion,3,2,4}(ncell,nid);
                CosSim{nregion,kk}(ncell,nid)=dot(A,B)/(norm(A)*norm(B));
                InnerP{nregion,kk}(ncell,nid)=dot(A,B);
            end
        end
        Ave(nregion,kk,1:4)=nanmean(CosSim{nregion,kk});
    end
    PermutationAve(:,:,:,nshuffle)=Ave;
end
save([OutPath,'/ShufflCosSim.mat'],'PermutationAve')

load([OutPath,'/3Period_CosSimilarity.mat'],'Ave')
for nregion=1:3
    P=squeeze(PermutationAve(nregion,3,:,:));
    P=sum(P);
    P=P(:);

    real=nanmean(Ave(nregion,3,:));

    p(nregion)=(length(find(P>real))+1)/101;
end