function [EBundle,measure,S] = bundleBasedOnRef(ParSet,extractedEndmember)
% endmember bundle extraction evaluation method
% the cluster of the bundle is based on the ref(reference), ref is a cell
% array that each element is one class of endmembers
% the EBundle is a cell array, each dimension has the same class of
% endmembers with the corresonding dimension of ref
% the result is a cell array that contains all evaluation results

img2d = ParSet.img2d;
ref = ParSet.endmemberBundle;
P = size(ref,1);
[L,numE] = size(extractedEndmember);

EBundle = cell(P,1);
bundleIndex = zeros(numE,1);
sad = zeros(P,1);
adm = zeros(P,1);
ads = zeros(P,1);
numEachClass = zeros(P,1);

%%  Find which class each endmember belongs to based on the set of reference endmembers
sadEachEndmember = zeros(numE,1);
for i = 1:numE
    tempMinSadEachClass = zeros(P,1);
    for j = 1:P
        tempRef = ref{j,1};
        numRefEachClass = size(tempRef,2);
        tempSadOneClass = zeros(numRefEachClass,1);
        for k = 1:numRefEachClass
            tempSadOneClass(k)= real(SAM(extractedEndmember(:,i),tempRef(:,k)));
        end
        tempMinSadEachClass(j) = min(tempSadOneClass);
    end
    [sadEachEndmember(i),bundleIndex(i)] = min(tempMinSadEachClass);
end

%% calculate the mean sad value for each class
for i = 1:P
    idx = bundleIndex == i;
    sad(i) = mean(sadEachEndmember(idx),'omitnan');
end

%% calculate the number of endmembers for each class
for i = 1:P
    idx = find(bundleIndex == i);
    numEachClass(i) = length(idx);
    EBundle{i} = extractedEndmember(:,idx);
end

%% calculate the mean difference of the mean spectral value from each class
 % of the estimated bundles and ref
endmemberForCalAbundance = zeros(L,P);
for i = 1:P
    meanEBundle = mean(EBundle{i},2,'omitnan');
    meanRef = mean(ref{i},2);
    adm(i) = sum(abs(meanEBundle-meanRef))/L;
    
    stdEBundle = std(EBundle{i},1,2,'omitnan'); 
    stdRef = std(ref{i},1,2);
    ads(i) = sum(abs(stdEBundle-stdRef))/L;
    
    endmemberForCalAbundance(:,i) = meanEBundle;
end

%% Calculating the RMSE
if ParSet.RMSEtf
    rmse = 100*ones(P,1);
    SBundle = sunsal([extractedEndmember,0.01*ones(L,1)],img2d,'POSITIVITY','yes','VERBOSE','yes','ADDONE','no', ...
        'LAMBDA', ParSet.LAMBDA,'AL_ITERS',ParSet.AL_ITERS);
    SBundle = SBundle(1:numE,:)./repmat(sum(SBundle(1:numE,:))+eps,numE,1);   

    S = zeros(P,ParSet.N);
    for i = 1:P
        idx = bundleIndex == i;
        S(i,:) = sum(SBundle(idx,:),1);
    end
    for i=1:P
        rmse(i) = sqrt(sum((ParSet.abundance(i,:)- S(i,:)).^2)/ParSet.N);
    end
    %% Calculating the RE
    ImgRmse =  F_rmse(img2d, extractedEndmember,SBundle);

%% If RMSE and RE is not assessed
else
    for i=1:P
        rmse(i) = 100;
        ImgRmse = 100;
    end
    S = 0;
end

measure.bundleIndex = bundleIndex;
measure.sad = sad;
measure.rmse = rmse;
measure.ImgRmse = ImgRmse;
measure.adm = adm;
measure.ads = ads;
measure.numEachClass = numEachClass; 
end

