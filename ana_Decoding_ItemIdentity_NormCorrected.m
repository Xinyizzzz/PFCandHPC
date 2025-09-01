%% mimic paper The Geometry of Abstraction in the Hippocampus
% and Prefrontal Cortex

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
%% only for IP period, only consider 8 conditions (8 I-Cue)
%% decode item identity within each co-location

clear all
clc
GroupName = {'HPC','PRC','PHC','TE','lPFC','mPFC','alPFC','FEF','dPS','vPS','dmPFC','vmPFC'};
repetitionN = 5;
CellScreen = {'','_3','_correct','_correct_3','_IP_4','_IP_4_nonsig'};
Group = zeros(8,2);
CLocationt = [1;1;2;2;3;3;4;4];
Icuet = [1;2;1;2;1;2;1;2];
Group(:,1) = CLocationt;
Group(:,2) = Icuet;
Predict = cell(2,2);
PredictFit=cell(4,repetitionN,2);
for n = 1:2
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
        Predict{n,2} ='A';
    elseif Group(n,2) == 2
        Predict{n,2} ='B';
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
end

for nCellScreen = 5 %%specify cell screen
    %% only correct trials or not
    if nCellScreen <3
        NCorrect=0;
    elseif nCellScreen>2
        NCorrect=1;
    end

    for nregion = [1 5 6]
        CellList=importdata(['E:\PFC and MTL\cell name\',GroupName{nregion},'_decoding',CellScreen{nCellScreen},'.txt']);
        cell_number = length(CellList);

        result_path = ['E:\PFC and MTL\Results\Review\decoding\ItemIdentity\IP_Bar\',GroupName{nregion},'_all_5\',CellScreen{nCellScreen},'\'];
        mkdir(result_path);
        performanceREd = zeros(4,20);
        performanceAved = zeros(4,20);
        performanceRER = zeros(4,20);
        performanceAveR = zeros(4,20);
        TimeMatrix = zeros(2,20);
        % to prepare for the re-sampled test set
        conditionN = [4;2];
        AllCondition = zeros(4,2);
        AllCondition(:,1) = [1;2;3;4];
        AllCondition(:,2) = [1;2;1;2];
        TestCondition = cell(2,2);
        TestCondition(1:2,1:2) = Predict(1:2,1:2);
        BIN = [0,500];
        for nbin = 1
            Accuracy_REd = zeros(4,repetitionN);
            Accuracy_aved = zeros(4,repetitionN);
            Accuracy_RER = zeros(4,repetitionN);
            Accuracy_aveR = zeros(4,repetitionN);
            for nCL = 1:4
                for nrep = 1:repetitionN
                    %% first, partition data
                    Train = zeros(cell_number,8,5);
                    Test = zeros(cell_number,8,5);
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


                        %% only decode item identity within each co-location
                        Icue1 = spike_matrix_correct(:,4);
                        CLocation1 = mod(Icue1,4);
                        CLocation1(find(CLocation1==0),:)=4;
                        MM = find(CLocation1==nCL);
                        spike_matrix_correct = spike_matrix_correct(MM,:);
                        Icue = Icue1(MM,1);
                        CLocation = CLocation1(MM,1);
                        for n = 1:length(spike_matrix_correct(:,1))
                            if Icue(n)<5
                                Icue(n)=1;
                            elseif Icue(n)>4
                                Icue(n)=2;
                            end
                        end

                        ip = zeros(length(spike_matrix_correct(:,1)),1);
                        for ntrial = 1:length(spike_matrix_correct(:,1))
                            ICueIndex = 3500;
                            START = ICueIndex+100;
                            END =  START+900;
                            ip(ntrial,1) = mean(spike_matrix_correct(ntrial,START:END),2)*1000;
                        end

                        k=1;
                        p=1;
                        for m = 1:2
                            A = find(Icue==Group(m,2));
                            Temp = randperm(length(A));
                            A = A(Temp',:);
                            Middle = floor(length(A)/2);
                            trainIndex = A(1:Middle,:);
                            testIndex = A((Middle+1):end,:);
                            Train(ncell,k:k+length(trainIndex)-1,1) = ip(trainIndex,:);
                            Train(ncell,k:k+length(trainIndex)-1,2) = CLocation(trainIndex,:);
                            Train(ncell,k:k+length(trainIndex)-1,3) = Icue(trainIndex,:);
                            Test(ncell,p:p+length(testIndex)-1,1)=ip(testIndex,:);
                            Test(ncell,p:p+length(testIndex)-1,2)=CLocation(testIndex,:);
                            Test(ncell,p:p+length(testIndex)-1,3)=Icue(testIndex,:);
                            k = k+length(trainIndex);
                            p = p+length(trainIndex);
                        end

                        % normalized by training set
                        AVE=mean(Train(ncell,:,1));
                        STD=std(Train(ncell,:,1));
                        Train(ncell,:,1)=(Train(ncell,:,1)-AVE)./STD;
                        Test(ncell,:,1)=(Test(ncell,:,1)-AVE)./STD;
                    end


                    %% second, re-sample
                    % re-sample 100 times for each condition for each trial
                    TrainRE = zeros(200,cell_number);
                    for nRE = 1:100
                        for m = 1:2
                            for ncell = 1:cell_number
                                AA = find(Train(ncell,:,3)==Group(m,2));
                                temp = randi(length(AA));
                                TrainRE((nRE-1)*2+m,ncell) = Train(ncell,AA(temp),1);
                            end
                        end
                    end
                    %         prepare inputTable
                    PredictRE = repmat(Predict,100,1);
                    TrainREsum = std(TrainRE);
                    B = find(TrainREsum==0|isnan(TrainREsum));
                    TrainRE(:,B)=[];
                    TrainRERaw = TrainRE;


                    %% third, test
                    % prepare re-sample
                    for ndecode = 2
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

                        %% form finalize
                        TrainRETemp = TrainRERaw;
                        TrainRETemp(:,C)=[];
                        TrainRE =  TrainRETemp;

                        TrainREs = num2cell(TrainRE);
                        TESTi = [TrainREs,PredictRE];
                        InputTablei_RE = cell2table(TESTi);

                        TestRE = num2cell(TestRE);
                        TESTi = TestRE;
                        testTableRE = cell2table(TESTi);

                        %% train and test
                        [trainedClassifier_RE, validationAccuracy] = LinearDiscriminant_ItemIdentity(InputTablei_RE,ndecode);
                        yfit_RER = trainedClassifier_RE.predictFcn(testTableRE);
                        correct_RER = 0;
                        for k = 1:size(yfit_RER,1)
                            if strcmp(yfit_RER{k,1},TestPredict{k,1})
                                correct_RER = correct_RER+1;
                            end
                        end
                        PredictFit{nCL,nrep,1}=yfit_RER;
                        PredictFit{nCL,nrep,2}=TestPredict;
                        Accuracy_RER(nCL,nrep) = correct_RER./size(yfit_RER,1);
                    end
                    disp([GroupName{nregion},'-',num2str(nbin),'-',num2str(nrep)])
                end
                save([result_path,'\Accuracy_RER_',num2str(nbin),'.mat'],'Accuracy_RER');
            end
            performanceRER(:,nbin) = mean(Accuracy_RER,2);
            %% write down each time bin
            TimeMatrix(1,nbin) = START;
            TimeMatrix(2,nbin) = END;
        end
        save([result_path,'\',GroupName{nregion},'_all.mat'],'performanceRER','PredictFit');
        save([result_path,'\TimeMatrix.mat'],'TimeMatrix');

    end
end