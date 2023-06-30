clc
%% Please cite our papers if you find it useful for your research.
% @ARTICLE{10158362,
%   author={Liu, Rong and Wang, Pengrui and Du, Bo and Qu, Bing},
%   journal={IEEE Geoscience and Remote Sensing Letters}, 
%   title={Endmember Bundle Extraction Based on Improved Multi-objective Particle Swarm Optimization}, 
%   year={2023},
%   volume={},
%   number={},
%   pages={1-1},
%   doi={10.1109/LGRS.2023.3287919}}

%% Read me
% 1.If you would like to test this code, 
% please refer to parset.m for how the data is entered and how the parameters are set,
% and we will give two datasets for a demo.

% 2.Note the parameter ParSet.th, 
% if you set it to 0 then in the code run you will need to enter the appropriate 
% threshold in matlab's command window based on the displayed histogram.

clc,clear,close all
%% Initialization parameters
problem = parset();
ParSet = problem.Berlin();   % Initialize parameter set ParSet (adjust to replace data set)
% ParSet = problem.urban();

numRun = ParSet.numRun;
P = ParSet.P;
Pt = ParSet.Ptrue;
sad = zeros(numRun,Pt);
rmse = zeros(numRun,Pt);
adm = zeros(numRun,Pt);
ads = zeros(numRun,Pt);
usedTime = zeros(numRun,1);
num = zeros(numRun,Pt);

%% Algorithm test 
for i = 1:numRun
    tic;
    % Extraction of candidate endmembers
    [indicies] = CCEBE(ParSet);
    fprintf('Complete EBE\n');
    % Using spatial post-processing in SSEBE
    w = ParSet.zoneL;
    percent = ParSet.percent;
    th = ParSet.th;
    [indiciesRemovedMixedPixels] = SSEBE(ParSet.img3d, w, percent, indicies, th);
    extractedEndmember = ParSet.img2d(:,indiciesRemovedMixedPixels);
    usedTime(i) = toc;
    fprintf('End of %d run, start of accuracy evaluation\n',i);

    %% Accuracy evaluation
    [EBundle,measure,SS] = bundleBasedOnRef(ParSet,extractedEndmember);
    sad(i,:)  =  measure.sad;
    rmse(i,:) =  measure.rmse;
    adm(i,:)  =  measure.adm;
    ads(i,:)  =  measure.ads;
    num(i,:)  =  measure.numEachClass;
    ImgRmse = measure.ImgRmse;
end

%% Output results
resultEachClass = {'meanSad','stdSad','meanRmse','stdRmse','meanAdm','stdAdm','meanAds','stdAds',...,
    'meanTime','stdTime','meanNum','stdNum','ImgRmse';mean(sad,'omitnan')',std(sad,'omitnan')',...,
    mean(rmse,'omitnan')',std(rmse,'omitnan')',mean(adm,'omitnan')',std(adm,'omitnan')',...,
    mean(ads,'omitnan')',std(ads,'omitnan')',mean(usedTime,'omitnan')',std(usedTime,'omitnan')',...,
    mean(num,'omitnan')',std(num,'omitnan')',ImgRmse};

fprintf('End of run\n');