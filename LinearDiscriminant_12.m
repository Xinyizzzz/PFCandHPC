function [trainedClassifier, validationAccuracy] = LinearDiscriminant(trainingData,decodedItem)
% [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% 返回经过训练的分类器及其准确度。以下代码重新创建在 Classification Learner App 中训
% 练的分类模型。您可以使用该生成的代码基于新数据自动训练同一模型，或通过它了解如何以程序化方
% 式训练模型。
%
%  输入:
%      trainingData: 一个所含预测变量和响应列与导入 App 中的相同的表。
%
%  输出:
%      trainedClassifier: 一个包含训练的分类器的结构体。该结构体中具有各种关于所训练分
%       类器的信息的字段。
%
%      trainedClassifier.predictFcn: 一个对新数据进行预测的函数。
%
%      validationAccuracy: 一个包含准确度百分比的双精度值。在 App 中，"历史记录" 列
%       表显示每个模型的此总体准确度分数。
%
% 使用该代码基于新数据来训练模型。要重新训练分类器，请使用原始数据或新数据作为输入参数
% trainingData 从命令行调用该函数。
%
% 例如，要重新训练基于原始数据集 T 训练的分类器，请输入:
%   [trainedClassifier, validationAccuracy] = trainClassifier(T)
%
% 要使用返回的 "trainedClassifier" 对新数据 T2 进行预测，请使用
%   yfit = trainedClassifier.predictFcn(T2)
%
% T2 必须是一个表，其中至少包含与训练期间使用的预测变量列相同的预测变量列。有关详细信息，请
% 输入:
%   trainedClassifier.HowToPredict

% 由 MATLAB 于 2022-05-06 22:10:51 自动生成


% 提取预测变量和响应
% 以下代码将数据处理为合适的形状以训练模型。
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

% 训练分类器
% 以下代码指定所有分类器选项并训练分类器。
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

% 使用预测函数创建结果结构体
predictorExtractionFcn = @(t) t(:, predictorNames);
discriminantPredictFcn = @(x) predict(classificationDiscriminant, x);
trainedClassifier.predictFcn = @(x) discriminantPredictFcn(predictorExtractionFcn(x));

% 向结果结构体中添加字段
trainedClassifier.RequiredVariables = predictorNames; %{'TESTi1', 'TESTi10', 'TESTi100', 'TESTi101', 'TESTi102', 'TESTi103', 'TESTi104', 'TESTi105', 'TESTi106', 'TESTi107', 'TESTi108', 'TESTi109', 'TESTi11', 'TESTi110', 'TESTi111', 'TESTi112', 'TESTi113', 'TESTi114', 'TESTi115', 'TESTi116', 'TESTi117', 'TESTi118', 'TESTi119', 'TESTi12', 'TESTi120', 'TESTi121', 'TESTi122', 'TESTi123', 'TESTi124', 'TESTi125', 'TESTi126', 'TESTi127', 'TESTi128', 'TESTi129', 'TESTi13', 'TESTi130', 'TESTi131', 'TESTi132', 'TESTi133', 'TESTi134', 'TESTi135', 'TESTi136', 'TESTi137', 'TESTi138', 'TESTi139', 'TESTi14', 'TESTi140', 'TESTi141', 'TESTi142', 'TESTi143', 'TESTi144', 'TESTi145', 'TESTi146', 'TESTi147', 'TESTi148', 'TESTi149', 'TESTi15', 'TESTi150', 'TESTi151', 'TESTi152', 'TESTi16', 'TESTi17', 'TESTi18', 'TESTi19', 'TESTi2', 'TESTi20', 'TESTi21', 'TESTi22', 'TESTi23', 'TESTi24', 'TESTi25', 'TESTi26', 'TESTi27', 'TESTi28', 'TESTi29', 'TESTi3', 'TESTi30', 'TESTi31', 'TESTi32', 'TESTi33', 'TESTi34', 'TESTi35', 'TESTi36', 'TESTi37', 'TESTi38', 'TESTi39', 'TESTi4', 'TESTi40', 'TESTi41', 'TESTi42', 'TESTi43', 'TESTi44', 'TESTi45', 'TESTi46', 'TESTi47', 'TESTi48', 'TESTi49', 'TESTi5', 'TESTi50', 'TESTi51', 'TESTi52', 'TESTi53', 'TESTi54', 'TESTi55', 'TESTi56', 'TESTi57', 'TESTi58', 'TESTi59', 'TESTi6', 'TESTi60', 'TESTi61', 'TESTi62', 'TESTi63', 'TESTi64', 'TESTi65', 'TESTi66', 'TESTi67', 'TESTi68', 'TESTi69', 'TESTi7', 'TESTi70', 'TESTi71', 'TESTi72', 'TESTi73', 'TESTi74', 'TESTi75', 'TESTi76', 'TESTi77', 'TESTi78', 'TESTi79', 'TESTi8', 'TESTi80', 'TESTi81', 'TESTi82', 'TESTi83', 'TESTi84', 'TESTi85', 'TESTi86', 'TESTi87', 'TESTi88', 'TESTi89', 'TESTi9', 'TESTi90', 'TESTi91', 'TESTi92', 'TESTi93', 'TESTi94', 'TESTi95', 'TESTi96', 'TESTi97', 'TESTi98', 'TESTi99'};
trainedClassifier.ClassificationDiscriminant = classificationDiscriminant;
trainedClassifier.About = '此结构体是从 Classification Learner R2019b 导出的训练模型。';
trainedClassifier.HowToPredict = sprintf('要对新表 T 进行预测，请使用: \n yfit = c.predictFcn(T) \n将 ''c'' 替换为作为此结构体的变量的名称，例如 ''trainedModel''。\n \n表 T 必须包含由以下内容返回的变量: \n c.RequiredVariables \n变量格式(例如矩阵/向量、数据类型)必须与原始训练数据匹配。\n忽略其他变量。\n \n有关详细信息，请参阅 <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>。');

% 提取预测变量和响应
% 以下代码将数据处理为合适的形状以训练模型。
%
inputTable = trainingData;
% predictorNames = %{'TESTi1', 'TESTi2', 'TESTi3', 'TESTi4', 'TESTi5', 'TESTi6', 'TESTi7', 'TESTi8', 'TESTi9', 'TESTi10', 'TESTi11', 'TESTi12', 'TESTi13', 'TESTi14', 'TESTi15', 'TESTi16', 'TESTi17', 'TESTi18', 'TESTi19', 'TESTi20', 'TESTi21', 'TESTi22', 'TESTi23', 'TESTi24', 'TESTi25', 'TESTi26', 'TESTi27', 'TESTi28'};
predictors = inputTable(:, predictorNames);
% response = inputTable.TESTi29;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];


% 计算重新代入的预测
[validationPredictions, validationScores] = predict(trainedClassifier.ClassificationDiscriminant, predictors);

% 计算重新代入的准确度
validationAccuracy = 1 - resubLoss(trainedClassifier.ClassificationDiscriminant, 'LossFun', 'ClassifError');
