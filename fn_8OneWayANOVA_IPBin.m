% function  Perc=fn_8OneWayANOVA_IPBin(InPath,OutPath)
InPath=['E:\PFC and MTL\Results\Review\DataMatrix2\'];
OutPath=['E:\PFC and MTL\Results\Review\ANOVA2\'];
load([InPath,'/DataMatrix.mat'],'DataMatrix');

Perc=nan(3,2);
P=cell(3,1);
FR=cell(3,1);

for nregion=1:size(DataMatrix,1)
    
    k1=0;k2=0;
    IPBin=DataMatrix(nregion).IPBin;
    Params=DataMatrix(nregion).Params;
    P{nregion}=nan(size(IPBin,2),2);
     FR{nregion}=nan(size(IPBin,2),1);
    
    for ncell=1:size(IPBin,2)
        Correct=Params{ncell}(:,end);
        Icue=Params{ncell}(:,1);
        indx=Correct>999&Icue>0;
        FR{nregion}(ncell,1)=DataMatrix(nregion).meanFR(ncell);
        
        % one-way ANOVA test
        p1=zeros(8,1);
        for nbin=1:8
            p1(nbin) =  anovan(IPBin{ncell}(indx,nbin), Params{ncell}(indx,1),'display','off');
        end
        P{nregion}(ncell,1)=min(p1);
        if min(p1)<0.00125&&FR{nregion}(ncell,1)>0.1
            k1=k1+1;
        end
        
        % kw test
        p2=zeros(8,1);
        for nbin=1:8
            p2(nbin) = kruskalwallis(IPBin{ncell}(indx,nbin), Params{ncell}(indx,1),'off');
        end
        P{nregion}(ncell,2)=min(p2);
        if min(p2)<0.00125&&FR{nregion}(ncell,1)>0.1
            k2=k2+1;
        end
    end
    
    Perc(nregion,1)=k1/size(IPBin,2);
    Perc(nregion,2)=k2/size(IPBin,2);
end

save([OutPath,'\8OneWay_IPBin.mat'],'Perc','P','FR')