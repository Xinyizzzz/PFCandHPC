function [trainedClassifier, validationAccuracy] = LinearDiscriminant(trainingData,decodedItem)
% [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% ���ؾ���ѵ���ķ���������׼ȷ�ȡ����´������´����� Classification Learner App ��ѵ
% ���ķ���ģ�͡�������ʹ�ø����ɵĴ�������������Զ�ѵ��ͬһģ�ͣ���ͨ�����˽�����Գ��򻯷�
% ʽѵ��ģ�͡�
%
%  ����:
%      trainingData: һ������Ԥ���������Ӧ���뵼�� App �е���ͬ�ı�
%
%  ���:
%      trainedClassifier: һ������ѵ���ķ������Ľṹ�塣�ýṹ���о��и��ֹ�����ѵ����
%       ��������Ϣ���ֶΡ�
%
%      trainedClassifier.predictFcn: һ���������ݽ���Ԥ��ĺ�����
%
%      validationAccuracy: һ������׼ȷ�Ȱٷֱȵ�˫����ֵ���� App �У�"��ʷ��¼" ��
%       ����ʾÿ��ģ�͵Ĵ�����׼ȷ�ȷ�����
%
% ʹ�øô��������������ѵ��ģ�͡�Ҫ����ѵ������������ʹ��ԭʼ���ݻ���������Ϊ�������
% trainingData �������е��øú�����
%
% ���磬Ҫ����ѵ������ԭʼ���ݼ� T ѵ���ķ�������������:
%   [trainedClassifier, validationAccuracy] = trainClassifier(T)
%
% Ҫʹ�÷��ص� "trainedClassifier" �������� T2 ����Ԥ�⣬��ʹ��
%   yfit = trainedClassifier.predictFcn(T2)
%
% T2 ������һ�����������ٰ�����ѵ���ڼ�ʹ�õ�Ԥ���������ͬ��Ԥ������С��й���ϸ��Ϣ����
% ����:
%   trainedClassifier.HowToPredict

% �� MATLAB �� 2022-05-06 22:10:51 �Զ�����


% ��ȡԤ���������Ӧ
% ���´��뽫���ݴ���Ϊ���ʵ���״��ѵ��ģ�͡�
%
inputTable = trainingData;
%  predictorNames = {'TESTi1', 'TESTi2', 'TESTi3', 'TESTi4', 'TESTi5', 'TESTi6', 'TESTi7', 'TESTi8', 'TESTi9', 'TESTi10', 'TESTi11', 'TESTi12', 'TESTi13', 'TESTi14', 'TESTi15', 'TESTi16', 'TESTi17', 'TESTi18', 'TESTi19', 'TESTi20', 'TESTi21', 'TESTi22', 'TESTi23', 'TESTi24', 'TESTi25', 'TESTi26', 'TESTi27', 'TESTi28', 'TESTi29', 'TESTi30', 'TESTi31', 'TESTi32', 'TESTi33', 'TESTi34', 'TESTi35', 'TESTi36', 'TESTi37', 'TESTi38', 'TESTi39', 'TESTi40', 'TESTi41', 'TESTi42', 'TESTi43', 'TESTi44', 'TESTi45', 'TESTi46', 'TESTi47', 'TESTi48', 'TESTi49', 'TESTi50', 'TESTi51', 'TESTi52', 'TESTi53', 'TESTi54', 'TESTi55', 'TESTi56', 'TESTi57', 'TESTi58', 'TESTi59', 'TESTi60', 'TESTi61', 'TESTi62', 'TESTi63', 'TESTi64', 'TESTi65', 'TESTi66', 'TESTi67', 'TESTi68', 'TESTi69', 'TESTi70', 'TESTi71', 'TESTi72', 'TESTi73', 'TESTi74', 'TESTi75', 'TESTi76', 'TESTi77', 'TESTi78', 'TESTi79', 'TESTi80', 'TESTi81', 'TESTi82', 'TESTi83', 'TESTi84', 'TESTi85', 'TESTi86', 'TESTi87', 'TESTi88', 'TESTi89', 'TESTi90', 'TESTi91', 'TESTi92', 'TESTi93', 'TESTi94', 'TESTi95', 'TESTi96', 'TESTi97', 'TESTi98', 'TESTi99', 'TESTi100', 'TESTi101', 'TESTi102', 'TESTi103', 'TESTi104', 'TESTi105', 'TESTi106', 'TESTi107', 'TESTi108', 'TESTi109', 'TESTi110', 'TESTi111', 'TESTi112', 'TESTi113', 'TESTi114', 'TESTi115', 'TESTi116', 'TESTi117', 'TESTi118', 'TESTi119', 'TESTi120', 'TESTi121', 'TESTi122', 'TESTi123', 'TESTi124', 'TESTi125', 'TESTi126', 'TESTi127', 'TESTi128', 'TESTi129', 'TESTi130', 'TESTi131', 'TESTi132', 'TESTi133', 'TESTi134', 'TESTi135', 'TESTi136', 'TESTi137', 'TESTi138', 'TESTi139', 'TESTi140', 'TESTi141', 'TESTi142', 'TESTi143', 'TESTi144', 'TESTi145', 'TESTi146', 'TESTi147', 'TESTi148', 'TESTi149', 'TESTi150', 'TESTi151', 'TESTi152'};
[nObservation,nTEST] = size(inputTable);
predictorNames = cell(1,nTEST-4);
for m = 1:nTEST-4
    predictorNames{1,m} = strcat('TESTi',num2str(m));
end
predictors = inputTable(:, predictorNames);
% response = inputTable.TESTi156;
response = inputTable(:,strcat('TESTi',num2str(nTEST-4+decodedItem)));
response = table2cell(response);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% ѵ��������
% ���´���ָ�����з�����ѡ�ѵ����������
if decodedItem == 1
    ClassNames = {'I';'II';'III';'IV'};
elseif decodedItem == 2
    ClassNames = {'1','2','3','4','5','6','7','8','9','10','11','12'};
elseif decodedItem==3
    ClassNames = {'-90','0','90'};
elseif decodedItem ==4
    ClassNames = {'TL';'TR';'BL';'BR'};
end
classificationDiscriminant = fitcdiscr(...
    predictors, ...
    response, ...
    'DiscrimType', 'diagLinear', ...
    'Gamma', 0, ...
    'FillCoeffs', 'off', ...
    'ClassNames', ClassNames);

% ʹ��Ԥ�⺯����������ṹ��
predictorExtractionFcn = @(t) t(:, predictorNames);
discriminantPredictFcn = @(x) predict(classificationDiscriminant, x);
trainedClassifier.predictFcn = @(x) discriminantPredictFcn(predictorExtractionFcn(x));

% �����ṹ��������ֶ�
trainedClassifier.RequiredVariables = predictorNames; %{'TESTi1', 'TESTi10', 'TESTi100', 'TESTi101', 'TESTi102', 'TESTi103', 'TESTi104', 'TESTi105', 'TESTi106', 'TESTi107', 'TESTi108', 'TESTi109', 'TESTi11', 'TESTi110', 'TESTi111', 'TESTi112', 'TESTi113', 'TESTi114', 'TESTi115', 'TESTi116', 'TESTi117', 'TESTi118', 'TESTi119', 'TESTi12', 'TESTi120', 'TESTi121', 'TESTi122', 'TESTi123', 'TESTi124', 'TESTi125', 'TESTi126', 'TESTi127', 'TESTi128', 'TESTi129', 'TESTi13', 'TESTi130', 'TESTi131', 'TESTi132', 'TESTi133', 'TESTi134', 'TESTi135', 'TESTi136', 'TESTi137', 'TESTi138', 'TESTi139', 'TESTi14', 'TESTi140', 'TESTi141', 'TESTi142', 'TESTi143', 'TESTi144', 'TESTi145', 'TESTi146', 'TESTi147', 'TESTi148', 'TESTi149', 'TESTi15', 'TESTi150', 'TESTi151', 'TESTi152', 'TESTi16', 'TESTi17', 'TESTi18', 'TESTi19', 'TESTi2', 'TESTi20', 'TESTi21', 'TESTi22', 'TESTi23', 'TESTi24', 'TESTi25', 'TESTi26', 'TESTi27', 'TESTi28', 'TESTi29', 'TESTi3', 'TESTi30', 'TESTi31', 'TESTi32', 'TESTi33', 'TESTi34', 'TESTi35', 'TESTi36', 'TESTi37', 'TESTi38', 'TESTi39', 'TESTi4', 'TESTi40', 'TESTi41', 'TESTi42', 'TESTi43', 'TESTi44', 'TESTi45', 'TESTi46', 'TESTi47', 'TESTi48', 'TESTi49', 'TESTi5', 'TESTi50', 'TESTi51', 'TESTi52', 'TESTi53', 'TESTi54', 'TESTi55', 'TESTi56', 'TESTi57', 'TESTi58', 'TESTi59', 'TESTi6', 'TESTi60', 'TESTi61', 'TESTi62', 'TESTi63', 'TESTi64', 'TESTi65', 'TESTi66', 'TESTi67', 'TESTi68', 'TESTi69', 'TESTi7', 'TESTi70', 'TESTi71', 'TESTi72', 'TESTi73', 'TESTi74', 'TESTi75', 'TESTi76', 'TESTi77', 'TESTi78', 'TESTi79', 'TESTi8', 'TESTi80', 'TESTi81', 'TESTi82', 'TESTi83', 'TESTi84', 'TESTi85', 'TESTi86', 'TESTi87', 'TESTi88', 'TESTi89', 'TESTi9', 'TESTi90', 'TESTi91', 'TESTi92', 'TESTi93', 'TESTi94', 'TESTi95', 'TESTi96', 'TESTi97', 'TESTi98', 'TESTi99'};
trainedClassifier.ClassificationDiscriminant = classificationDiscriminant;
trainedClassifier.About = '�˽ṹ���Ǵ� Classification Learner R2019b ������ѵ��ģ�͡�';
trainedClassifier.HowToPredict = sprintf('Ҫ���±� T ����Ԥ�⣬��ʹ��: \n yfit = c.predictFcn(T) \n�� ''c'' �滻Ϊ��Ϊ�˽ṹ��ı��������ƣ����� ''trainedModel''��\n \n�� T ����������������ݷ��صı���: \n c.RequiredVariables \n������ʽ(�������/��������������)������ԭʼѵ������ƥ�䡣\n��������������\n \n�й���ϸ��Ϣ������� <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>��');

% ��ȡԤ���������Ӧ
% ���´��뽫���ݴ���Ϊ���ʵ���״��ѵ��ģ�͡�
%
inputTable = trainingData;
% predictorNames = %{'TESTi1', 'TESTi2', 'TESTi3', 'TESTi4', 'TESTi5', 'TESTi6', 'TESTi7', 'TESTi8', 'TESTi9', 'TESTi10', 'TESTi11', 'TESTi12', 'TESTi13', 'TESTi14', 'TESTi15', 'TESTi16', 'TESTi17', 'TESTi18', 'TESTi19', 'TESTi20', 'TESTi21', 'TESTi22', 'TESTi23', 'TESTi24', 'TESTi25', 'TESTi26', 'TESTi27', 'TESTi28'};
predictors = inputTable(:, predictorNames);
% response = inputTable.TESTi29;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];


% �������´����Ԥ��
[validationPredictions, validationScores] = predict(trainedClassifier.ClassificationDiscriminant, predictors);

% �������´����׼ȷ��
validationAccuracy = 1 - resubLoss(trainedClassifier.ClassificationDiscriminant, 'LossFun', 'ClassifError');
