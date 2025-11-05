function fn_GetSpikeProperteis_fromPLX_v2(PlxPath,OutPath)

PlxPath={['E:\monkey\Care\PFC\plx\done\'];['E:\monkey\Galeno\PFC\plx\done\']};
OutPath=['E:\PFC and MTL\Results\Review\SpikeMatrix\'];
mkdir(OutPath)
CellName=cell(length(PlxPath),1);
CellInfo=cell(length(PlxPath),1);

for nsec=1:length(PlxPath)
    fileList = dir( fullfile(PlxPath{nsec}, '*.plx') );
    FileNames = {fileList.name};

    ncell=0;

    for nfile=1:length(FileNames)
        PlxName=[PlxPath{nsec},FileNames{nfile}];
        try
            [~,names] = plx_chan_names(PlxName);% channel name
        catch
            continue
        end
        tscounts = plx_info(PlxName, 1);
        unitnum=length(find(tscounts>0))-1;

        for unit=1:unitnum
            ncell=ncell+1;

            [nSpks, ts] = plx_ts(PlxName, names, unit);


            % firing rate
            fr=nSpks/(max(ts)-min(ts));

            % burst rate
            ISI=diff(ts);
            br=length(find(ISI<0.02))/length(ISI);

            % spike waveform width
            fs = 44000;
            meanWidth = SpkWfWidth(fs, PlxName, names, unit);

            % cell name
            cellname=[FileNames{nfile}(1:end-4),'-1-',num2str(unit)];


            CellName{nsec}{ncell,1}=cellname;
            CellInfo{nsec}(ncell,1)=fr;
            CellInfo{nsec}(ncell,2)=br;
            CellInfo{nsec}(ncell,3)=meanWidth;
        end
    end
    save([OutPath,'/SpikeProperties_v2.mat'],'CellInfo','CellName')
end

load('E:\PFC and MTL\Results\Review\SpikeMatrix\SpikeProperties_v2.mat','CellInfo','CellName')
CellInfo1=CellInfo;
CellName1=CellName;
clear CellInfo CellName

load('E:\PFC and MTL\Results\Review\HPC\HPC\Spike Properties\SpikeProperties_HPC.mat','CellName','CellInfo')
CellName2=CellName;
CellInfo2=CellInfo(:,1:3);

clear CellName CellInfo

CellInfo=CellInfo1;
CellName=CellName1;
CellInfo{3}=CellInfo2;
CellName{3}=CellName2;
save([OutPath,'/SpikeProperties_v2_full.mat'],'CellInfo','CellName')




% load('E:\PFC and MTL\Results\Review\SpikeMatrix\SpikeProperties.mat')
% CellInfo1=CellInfo;
% CellName1=CellName;
% 
% load('E:\PFC and MTL\Results\Review\SpikeMatrix\SpikeProperties_2.mat')
% CellInfo2=CellInfo;
% CellName2=CellName;
% 
% CellInfo{1}=CellInfo1{1};
% CellInfo{2}=cat(1,CellInfo1{2},CellInfo2{2});
% 
% CellName{1}=CellName1{1};
% CellName{2}=cat(1,CellName1{2},CellName2{2});
