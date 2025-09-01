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
repetitionN = 100;
PredictFit=cell(repetitionN,16,2);
CellScreen = {'','_3','_correct','_correct_3','_IP_4','_IP_4_nonsig'};
Group = zeros(8,3);
CLocationt = [1;2;3;4];
Icuet = [1;2;3;4;5;6;7;8];
Angle1 = -90+zeros(8,1);
Angle2 = zeros(8,1);
Angle3 = 90+zeros(8,1);
%     Shiftt = [1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;3;3;3;3;3;3;3;3];
Group(:,1) = repmat(CLocationt,[2,1]);
Group(:,2) = repmat(Icuet,[1,1]);
% Group(:,3) = [Angle1;Angle2;Angle3];
Predict = cell(8,4);

for n = 1:8
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
end

%% first of first, consider all possible combinations for training and testing
DesignMatrix(:,:,1) = [1,5;2,6;3,7;4,8];DesignMatrix(:,:,2) = [1,5;2,6;3,7;8,4];DesignMatrix(:,:,3) = [1,5;2,6;7,3;4,8];DesignMatrix(:,:,4) = [1,5;2,6;7,3;8,4];
DesignMatrix(:,:,5) = [1,5;6,2;3,7;4,8];DesignMatrix(:,:,6) = [1,5;6,2;3,7;8,4];DesignMatrix(:,:,7) = [1,5;6,2;7,3;4,8];DesignMatrix(:,:,8) = [1,5;6,2;7,3;8,4];
DesignMatrix(:,:,9) = [5,1;2,6;3,7;4,8];DesignMatrix(:,:,10) = [5,1;2,6;3,7;8,4];DesignMatrix(:,:,11) =  [5,1;2,6;7,3;4,8];DesignMatrix(:,:,12) = [5,1;2,6;7,3;8,4];
DesignMatrix(:,:,13) = [5,1;6,2;3,7;4,8];DesignMatrix(:,:,14) = [5,1;6,2;3,7;4,8];DesignMatrix(:,:,15) = [5,1;6,2;3,7;4,8];DesignMatrix(:,:,16) = [5,1;6,2;3,7;4,8];


% to prepare for the re-sampled test set
conditionN = [4;8;3;4];
AllCondition = zeros(8,4);
AllCondition(:,1) = [1;2;3;4;1;2;3;4];
AllCondition(:,2) = [1;2;3;4;5;6;7;8];
AllCondition(:,3) = [-90;0;90;-90;0;90;0;0];
AllCondition(:,4) = [1;2;3;4;1;2;3;4];
TestCondition = cell(8,4);
TestCondition(1:8,1:2) = Predict(1:8,1:2);
TestCondition(:,3) = {'-90','0','90','-90','0','90','0','0'};
TestCondition(:,4) = {'TR','TL','BL','BR','TR','TL','BL','BR'};
for nCellScreen = 5 %%specify cell screen
    %% only correct trials or not
    if nCellScreen <3
        NCorrect=0;
    elseif nCellScreen>2
        NCorrect=1;
    end
    
    for nregion = [1 5 6]
        CellList=importdata(['E:\PFC and MTL\cell name\',GroupName{nregion},'_decoding',CellScreen{nCellScreen},'.txt']);
        %         CellList=importdata(['E:\PFC and MTL\cell name\cell name temp.txt']);
        cell_number = length(CellList);
        
        result_path = ['E:\PFC and MTL\Results\Review\decoding\ABSet_4\IP_Bar\',GroupName{nregion},'_all_AllComb_100\',CellScreen{nCellScreen},'\'];
        mkdir(result_path);
        performanceREd = zeros(4,20);
        performanceAved = zeros(4,20);
        performanceRER = zeros(16,20);
        performanceAveR = zeros(4,20);
        TimeMatrix = zeros(2,20);
        
        for nbin = 1
            Accuracy_REd = zeros(4,repetitionN);
            Accuracy_aved = zeros(4,repetitionN);
            Accuracy_RER = zeros(16,repetitionN);
            Accuracy_aveR = zeros(4,repetitionN);
            for nrep = 1:repetitionN
                
                for ncomb=1:16
                    
                    TrainGroup = Group(DesignMatrix(:,1,ncomb),:);
                    TrainPredict = Predict(DesignMatrix(:,1,ncomb),:);
                    TestGroup = Group(DesignMatrix(:,2,ncomb),:);
                    TestPredict = Predict(DesignMatrix(:,2,ncomb),:);
                    
                    
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
                        for ntrial = 1:length(spike_matrix_correct(:,1))
                            ICueIndex = 3500;
                            START = ICueIndex+100;
                            END =  START+900;
                            ip(ntrial,1) = mean(spike_matrix_correct(ntrial,START:END),2)*1000;
                        end
                        
                        
                        k=1;
                        p=1;
                        for m = 1:4
                            A = find(Icue==TrainGroup(m,2));
                            trainIndex = A;
                            Train(ncell,k:k+length(trainIndex)-1,1) = ip(trainIndex,:);
                            Train(ncell,k:k+length(trainIndex)-1,2) = CLocation(trainIndex,:);
                            Train(ncell,k:k+length(trainIndex)-1,3) = Icue(trainIndex,:);
                            Train(ncell,k:k+length(trainIndex)-1,4) = Ocue(trainIndex,:);
                            Train(ncell,k:k+length(trainIndex)-1,5) = Tar(trainIndex,:);
                            B = find(Icue==TestGroup(m,2));
                            testIndex = B;
                            Test(ncell,p:p+length(testIndex)-1,1)=ip(testIndex,:);
                            Test(ncell,p:p+length(testIndex)-1,2)=CLocation(testIndex,:);
                            Test(ncell,p:p+length(testIndex)-1,3)=Icue(testIndex,:);
                            Test(ncell,p:p+length(testIndex)-1,4)=Ocue(testIndex,:);
                            Test(ncell,p:p+length(testIndex)-1,5)=Tar(testIndex,:);
                            k = k+length(trainIndex);
                            p = p+length(testIndex);
                            clear A B
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
                    TrainRE = zeros(400,cell_number);
                    for nRE = 1:100
                        for m = 1:4
                            for ncell = 1:cell_number
                                AA = find(Train(ncell,:,3)==TrainGroup(m,2));
                                temp = randi(length(AA));
                                TrainRE((nRE-1)*4+m,ncell) = Train(ncell,AA(temp),1);
                            end
                        end
                    end
                    %         prepare inputTable
                    PredictRE = repmat(TrainPredict,100,1);
                    TrainREsum = std(TrainRE);
                    B = find(TrainREsum==0|isnan(TrainREsum));
                    TrainRE(:,B)=[];
                    TrainRERaw = TrainRE;
                    
                    
                    
                    %% third, test
                    % prepare re-sample
                    for ndecode = 1
                        TestRE = zeros(conditionN(ndecode,1)*100,cell_number);
                        TestPredictRE = cell(conditionN(ndecode,1)*100,1);
                        for nRE = 1:100
                            for ncon = 1:conditionN(ndecode,1)
                                for ncell = 1:cell_number
                                    AAA = find(Test(ncell,:,ndecode+1)==AllCondition(ncon,ndecode));
                                    temp = randi(length(AAA));
                                    TestRE((nRE-1)*conditionN(ndecode,1)+ncon,ncell) = Test(ncell,AAA(temp),1);
                                end
                                TestPredictRE{(nRE-1)*conditionN(ndecode,1)+ncon,1} = TestCondition{ncon,ndecode};
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
                        [trainedClassifier_RE, validationAccuracy] = LinearDiscriminant(InputTablei_RE,ndecode);
                        yfit_RER = trainedClassifier_RE.predictFcn(testTableRE);
                        correct_RER = 0;
                        for k = 1:size(yfit_RER,1)
                            if strcmp(yfit_RER{k,1},TestPredictRE{k,1})
                                correct_RER = correct_RER+1;
                            end
                        end
                        PredictFit{nrep,ncomb,1}=yfit_RER;
                        PredictFit{nrep,ncomb,2}=TestPredictRE;
                        Accuracy_RER(ncomb,nrep) = correct_RER./size(yfit_RER,1);
                        %test on re-sampled test set
                        % wait to be finished
                    end
                    disp([GroupName{nregion},'-',num2str(ncomb),'-',num2str(nrep)])
                end
            end
            save([result_path,'\Accuracy_RER_',num2str(nbin),'.mat'],'Accuracy_RER');
            performanceRER(:,nbin) = mean(Accuracy_RER,2);
            %% write down each time bin
            TimeMatrix(1,nbin) = START;
            TimeMatrix(2,nbin) = END;
            %     performanceAveR(:,nbin) = mean(Accuracy_aveR,2);
        end
        save([result_path,'\',GroupName{nregion},'_all.mat'],'performanceRER','PredictFit','-v7.3');
        save([result_path,'\TimeMatrix.mat'],'TimeMatrix');
        AAAAA(nregion,1) = performanceRER(1,1);
    end
end