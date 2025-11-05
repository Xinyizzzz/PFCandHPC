% function [Perc,PercShuffle]=fn_ThreeeWayNested_OPBin(InPath,OutPath)
 InPath=['E:\PFC and MTL\Results\Review\DataMatrix\'];
    OutPath=['E:\PFC and MTL\Results\Review\ANOVA\'];
load([InPath,'/DataMatrix.mat'],'DataMatrix');

Perc=nan(3,4);
PercShuffle=nan(3,4);
SigIdentity=cell(3,2);
FR=cell(3,1);

for nregion=1:size(DataMatrix,1)

    OPBin=DataMatrix(nregion).OPBin;
    Params=DataMatrix(nregion).Params;
    %Params=[Icue,CL,Ocue,Tar,Choice,Correct];
    Nsig=zeros(size(OPBin,2),4);
    Nsigshuffle=zeros(size(OPBin,2),4);
    FR{nregion}=nan(size(OPBin,2),1);
    

    for ncell=1:size(OPBin,2)
        disp([num2str(nregion),'-',num2str(ncell)])
        Correct=Params{ncell}(:,end);
        Tar=Params{ncell}(indx,4);
        Ocue=Params{ncell}(indx,3);
        Icue2=Params{ncell}(indx,2);
        Icue=Params{ncell}(indx,1);
        indx=Correct>999;
        F=nan(8,4);np=zeros(8,4);
        
        FR{nregion}(ncell,1)=DataMatrix(nregion).meanFR(ncell);
        
        if FR{nregion}(ncell,1)<=0.1
            continue
        end
        
        for i = 1:8
             Act=OPBin{ncell}(indx,i);
            [p,table] = anovan(Act,{Tar,Ocue,Icue2,Icue},'display','off','nested',[0,0,0,0;0,0,0,0;0,0,0,0;0,0,1,0],'varnames',{'target','ocue','co-location','icue'});
            np(i,:) = p';
            F(i,1)=table{2,6};
            F(i,2)=table{3,6};
            F(i,3)=table{4,6};
            F(i,4)=table{5,6};
        end
        % ANOVA
        pvalue(1,1) = min(np(:,1));
        pvalue(1,2) = min(np(:,2));
        pvalue(1,3) = min(np(:,3));
        pvalue(1,4) = min(np(:,4));
        for nvar = 1:4
            if pvalue(1,nvar) <0.00125
                Nsig(ncell,nvar) =1;
            end
        end
        OPBinNow=OPBin{ncell};
        % permutation
        Fshuffle=nan(8,4,2000);
        parfor nshuffle=1:2000 %1000 is the minimum to reach 0.00125
            Ftemp=nan(8,4);
            for i = 1:8
                Act=OPBinNow(indx,i);
                Actshuffle=Act(randperm(length(Act)));
                [~,table] = anovan(Actshuffle,{Tar,Ocue,Icue2,Icue},'display','off','nested',[0,0,0,0;0,0,0,0;0,0,0,0;0,0,1,0],'varnames',{'target','ocue','co-location','icue'});
                Ftemp(i,1)=table{2,6};
                Ftemp(i,2)=table{3,6};
                Ftemp(i,3)=table{4,6};
                Ftemp(i,4)=table{5,6};
            end
            Fshuffle(:,:,nshuffle)=Ftemp;
        end
        FAll{nregion}{ncell,1}=F;
        FAll{nregion}{ncell,2}=Fshuffle;
        pshuffle=nan(8,4);
        for nvar=1:4
            for nbin=1:8
                if ~isnan(F(nbin,nvar))
                pshuffle(nbin,nvar)=(length(find(Fshuffle(nbin,nvar,:)>F(nbin,nvar)))+1)/2001;
                else
                    pshuffle(nbin,nvar)=nan;
                end
            end
            if min(pshuffle(:,nvar))<=0.00125
                Nsigshuffle(ncell,nvar)=1;
            end
        end
    end
    
    Perc(nregion,1:4)=sum(Nsig)./size(OPBin,2);
    PercShuffle(nregion,1:4)=sum(Nsigshuffle)./size(OPBin,2);
    SigIdentity{nregion,1}=Nsig;
    SigIdentity{nregion,2}=Nsigshuffle;
end

save([OutPath,'\ThreeWayNested_wShuffle_OPBin_3.mat'],'Perc','PercShuffle','SigIdentity','FAll')

clear np p
