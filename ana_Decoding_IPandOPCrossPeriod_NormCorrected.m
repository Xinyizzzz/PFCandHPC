%% mimic paper The Geometry of Abstraction in the Hippocampus
% and Prefrontal Cortex]

%% for each repeatition
%% first
% for each contiditon and each neuron, leave one trial as test set, at
% least two trials as train set. Cells that doesn't have three repetitions
% in each condition should be excluded, use ana_CellScreen.
%% second
% for training set, normalize for each cell across all conditions.
% (1) average within each consition
% (2) 1000 re-sample for each condition for each neuron
%% third
% for test set, need to be normalized or not ?
% should be re-sampled, but not neccessary because only have one trial for
% each condition.(1) directly test
% another way is to re-sample 1000 times based on each decoded item.
% 1-co-location,2-Icue,3-ocue,4-target
clear all
clc
GroupName = {'HPC','PRC','PHC','TE','lPFC','mPFC','alPFC','FEF','dPS','vPS','dmPFC','vmPFC'};
repetitionN = 900;
PredictFit=cell(repetitionN,2);
CellScreen = {'','_3','_correct','_correct_3'};
Group = zeros(24,3);
CLocationt = [1;2;3;4];
Icuet = [1;2;3;4;5;6;7;8];
Angle1 = -90+zeros(8,1);
Angle2 = zeros(8,1);
Angle3 = 90+zeros(8,1);
%     Shiftt = [1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;3;3;3;3;3;3;3;3];
Group(:,1) = repmat(CLocationt,[6,1]);
Group(:,2) = repmat(Icuet,[3,1]);
Group(:,3) = [Angle1;Angle2;Angle3];
Predict = cell(24,4);
for n = 1:24
    %col-location
    if Group(n,1) == 1
        Predict{n,1} = 'I';
    elseif Group(n,1) == 2
        Predict{n,1} = 'II';
    elseif Group(n,1) == 3
        Predict{n,1} = 'III';
    elseif Group(n,1)==4
        Predict{n,1}='IV';
    end
    
    %icue
    if Group(n,2) == 1
        Predict{n,2} ='IA';
    elseif Group(n,2) == 2
        Predict{n,2} ='IIA';
    elseif Group(n,2) == 3
        Predict{n,2} ='IIIA';
    elseif Group(n,2) == 4
        Predict{n,2} ='IVA';
    elseif Group(n,2) == 5
        Predict{n,2} ='IB';
    elseif Group(n,2) == 6
        Predict{n,2} ='IIB';
    elseif Group(n,2) == 7
        Predict{n,2} ='IIIB';
    elseif Group(n,2) == 8
        Predict{n,2} ='IVB';
    end
    
    %Ocue
    Predict{n,3} = num2str(Group(n,3));
    
    % target
    if Group(n,1)==1
        if Group(n,3)==-90
            Predict{n,4} = 'TL';
        elseif Group(n,3)==0
            Predict{n,4} = 'TR';
        elseif Group(n,3)==90
            Predict{n,4} = 'BR';
        end
    elseif Group(n,1)==2
        if Group(n,3)==-90
            Predict{n,4} = 'BL';
        elseif Group(n,3)==0
            Predict{n,4} = 'TL';
        elseif Group(n,3)==90
            Predict{n,4} = 'TR';
        end
    elseif Group(n,1)==3
        if Group(n,3)==-90
            Predict{n,4} = 'BR';
        elseif Group(n,3)==0
            Predict{n,4} = 'BL';
        elseif Group(n,3)==90
            Predict{n,4} = 'TL';
        end
    elseif Group(n,1)==4
        if Group(n,3)==-90
            Predict{n,4} = 'TR';
        elseif Group(n,3)==0
            Predict{n,4} = 'BR';
        elseif Group(n,3)==90
            Predict{n,4} = 'BL';
        end
    end
end

for nCellScreen = 4 %%specify cell screen
    %% only correct trials or not
    if nCellScreen <3
        NCorrect=0;
    elseif nCellScreen>2
        NCorrect=1;
    end
    
    for nregion = [1 5 6]
        CellList=importdata(['E:\PFC and MTL\cell name\',GroupName{nregion},'_decoding',CellScreen{nCellScreen},'.txt']);
        cell_number = length(CellList);
        
        result_path = ['E:\PFC and MTL\Results\Review\decoding\IPandOPCrossPeriod\Bar_OP\',GroupName{nregion},'_all_900\',CellScreen{nCellScreen},'\'];
        mkdir(result_path);
        performanceREd = zeros(4,20);
        performanceAved = zeros(4,20);
        performanceRER = zeros(4,20);
        performanceAveR = zeros(4,20);
        TimeMatrix = zeros(2,20);
        
        % to prepare for the re-sampled test set
        conditionN = [4;8;3;4];
        AllCondition = zeros(8,4);
        AllCondition(:,1) = [1;2;3;4;1;2;3;4];
        AllCondition(:,2) = [1;2;3;4;5;6;7;8];
        AllCondition(:,3) = [-90;0;90;-90;0;90;0;0];
        if nregion < 5
            AllCondition(:,4) = [1;3;5;7;1;3;5;7;];
        elseif nregion > 4
            AllCondition(:,4) = [1;2;3;4;1;2;3;4];
        end
        TestCondition = cell(8,4);
        TestCondition(1:8,1:2) = Predict(1:8,1:2);
        TestCondition(:,3) = {'-90','0','90','-90','0','90','0','0'};
        TestCondition(:,4) = {'I','II','III','IV','I','II','III','IV'};
        TestCondition(:,4) = {'TR','TL','BL','BR','TR','TL','BL','BR'};
        
        for nbin1 = 1
            for nbin2 = 1
                Accuracy_REd = zeros(4,repetitionN);
                Accuracy_aved = zeros(4,repetitionN);
                Accuracy_RER = zeros(1,repetitionN);
                Accuracy_aveR = zeros(4,repetitionN);
                for nrep = 1:repetitionN
                    
                    %% first, partition data
                    Train = zeros(cell_number,24,5);
                    Test = zeros(cell_number,24,5);
                    for ncell = 1:cell_number
                        CellName = CellList{ncell,1};
                        global_path = 'E:\PFC and MTL\spike matrix\';
                        load([global_path,CellName,'_spike.mat']);
                        
                        %% correct trials only
                        if NCorrect==1
                            spike_matrix_all = spike_matrix_fin;
                            if nregion <5
                                spike_matrix_correct = spike_matrix_all(find(spike_matrix_all(:,7)>20000),:);
                            elseif nregion >4
                                spike_matrix_correct = spike_matrix_all(find(spike_matrix_all(:,7)>999&spike_matrix_all(:,4)>0),:);% only correct trials
                            end
                        elseif NCorrect==0
                            spike_matrix_correct = spike_matrix_fin;
                        end
                        
                        %                         %% Shuffle  within co-location
                        %                         spike_matrix_fin_raw = spike_matrix_fin;
                        %                         Tar = mod(spike_matrix_fin(:,9),10);
                        %                         spike_matrix_fin(:,35) = Tar;
                        %                         spike_matrix_fin = sortrows(spike_matrix_fin,35);
                        %                         for cl = 1:8
                        %                             A = find(spike_matrix_fin(:,35)==cl);
                        %                             M = spike_matrix_fin(A,:);
                        %                             Temp = randperm(length(M(:,1)));
                        %                             MM = M(Temp',:);
                        %                             spike_matrix_fin(A,1:20) = MM(:,1:20);
                        %                         end
                        
                        Icue = spike_matrix_correct(:,4);
                        CLocation = mod(Icue,4);
                        CLocation(find(CLocation==0),:)=4;
                        Ocue = spike_matrix_correct(:,5);
                        if nregion > 4
                            Tar = mod(spike_matrix_correct(:,9),10);
                        elseif nregion < 5
                            Tar = spike_matrix_correct(:,8);
                        end
                        
                        ip = zeros(length(spike_matrix_correct(:,1)),1);
                        op = zeros(length(spike_matrix_correct(:,1)),1);
                        for ntrial = 1:length(spike_matrix_correct(:,1))
                            if nregion < 5
                                OCueIndex = spike_matrix_fin(ntrial,16)+1000;
                            elseif nregion > 4
                                OCueIndex = spike_matrix_fin(ntrial,19);
                            end
                            ICueIndex = 3500;
                            START1 = ICueIndex+(nbin1-1)*300+100;
                            END1 =  START1+900;
                            START2 = OCueIndex+(nbin2-1)*300+100;
                            END2 =  START2+900;
                            ip(ntrial,1) = mean(spike_matrix_correct(ntrial,START1:END1),2)*1000;
                            op(ntrial,1) = mean(spike_matrix_correct(ntrial,START2:END2),2)*1000;
                        end
                        
                        
                        k=1;
                        for m = 1:24
                            A = find(Icue==Group(m,2)&Ocue==Group(m,3));
                            testIndex = randi(length(A));
                            trainIndex = A;
                            trainIndex(testIndex,:)=[];
                            testIndex = A(testIndex,:);
                            Train(ncell,k:k+length(trainIndex)-1,1) = ip(trainIndex,:);
                            Train(ncell,k:k+length(trainIndex)-1,2) = CLocation(trainIndex,:);
                            Train(ncell,k:k+length(trainIndex)-1,3) = Icue(trainIndex,:);
                            Train(ncell,k:k+length(trainIndex)-1,4) = Ocue(trainIndex,:);
                            Train(ncell,k:k+length(trainIndex)-1,5) = Tar(trainIndex,:);
                            Test(ncell,m,1)=op(testIndex,:);
                            Test(ncell,m,2)=CLocation(testIndex,:);
                            Test(ncell,m,3)=Icue(testIndex,:);
                            Test(ncell,m,4)=Ocue(testIndex,:);
                            Test(ncell,m,5)=Tar(testIndex,:);
                            k = k+length(trainIndex);
                        end
                        
                        % normalized by training set
                        AVE=mean(Train(ncell,:,1));
                        STD=std(Train(ncell,:,1));
                        Train(ncell,:,1)=(Train(ncell,:,1)-AVE)./STD;
                        % apply training set normalization parameters to
                        % testing set
                        Test(ncell,:,1)=(Test(ncell,:,1)-AVE)./STD;
                    end
                    
                    %% second, re-sample
                    % re-sample 100 times for each condition for each trial
                    TrainRE = zeros(2400,cell_number);
                    for nRE = 1:100
                        for m = 1:24
                            for ncell = 1:cell_number
                                AA = find(Train(ncell,:,3)==Group(m,2)&Train(ncell,:,4)==Group(m,3));
                                temp = randi(length(AA));
                                TrainRE((nRE-1)*24+m,ncell) = Train(ncell,AA(temp),1);
                            end
                        end
                    end
                    %         prepare inputTable
                    PredictRE = repmat(Predict,100,1);
                    TrainREsum = std(TrainRE);
                    B = find(TrainREsum==0|isnan(TrainREsum));
                    TrainRE(:,B)=[];
                    
                    
                    
                    %% third, test
                    % prepare re-sample
                    for ndecode = 1
                        TestRE = zeros(conditionN(ndecode,1)*100,cell_number);
                        TestPredict = cell(conditionN(ndecode,1)*100,1);
                        for nRE = 1:100
                            for ncon = 1:conditionN(ndecode,1)
                                for ncell = 1:cell_number
                                    AAA = find(Test(ncell,:,ndecode+1)==AllCondition(ncon,ndecode));
                                    temp = randi(length(AAA));
                                    TestRE((nRE-1)*conditionN(ndecode,1)+ncon,ncell) = Test(ncell,AAA(temp),1);
                                end
                                TestPredict{(nRE-1)*conditionN(ndecode,1)+ncon,1} = TestCondition{ncon,ndecode};
                            end
                        end
                        TestRE(:,B)=[];
                        TestREsum = std(TestRE);
                        C = find(TestREsum==0|isnan(TestREsum));
                        TestRE(:,C)=[];
                        
                        %% form transformation
                        TrainRE(:,C)=[];
                        TrainREs = num2cell(TrainRE);
                        TESTi = [TrainREs,PredictRE];
                        InputTablei_RE = cell2table(TESTi);
                        
                        TestRE = num2cell(TestRE);
                        TESTi = TestRE;
                        testTableRE = cell2table(TESTi);
                        
                        %% train and test
                        [trainedClassifier_RE, validationAccuracy] = LinearDiscriminant(InputTablei_RE,ndecode);
                        yfit_RER = trainedClassifier_RE.predictFcn(testTableRE);
                        correct_RER = 0;
                        for k = 1:size(yfit_RER,1)
                            if strcmp(yfit_RER{k,1},TestPredict{k,1})
                                correct_RER = correct_RER+1;
                            end
                        end
                        PredictFit{nrep,1}=yfit_RER;
                        PredictFit{nrep,2}=TestPredict;
                        Accuracy_RER(ndecode,nrep) = correct_RER./size(yfit_RER,1);
                    end
                    disp([GroupName{nregion},'-',num2str(nbin1),'-',num2str(nbin2),'-',num2str(nrep)])
                end
                save([result_path,'\Accuracy_RER_',num2str(nbin1),'-',num2str(nbin2),'.mat'],'Accuracy_RER');
                performanceRER(nbin1,nbin2) = mean(Accuracy_RER(ndecode,:));
                %% write down each time bin
                TimeMatrix(1,nbin1) = START1;
                TimeMatrix(2,nbin1) = END1;
            end
        end
        save([result_path,'\',GroupName{nregion},'_all.mat'],'performanceRER');
        save([result_path,'\TimeMatrix.mat'],'TimeMatrix');
        AAAAA(nregion,1) = performanceRER(1,1);
    end
end