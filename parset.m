%% You need to import the data here and set the parameters

function [problem] = parset()

    problem.Berlin = @Berlin;
    problem.urban = @urban;
    
end

%% BerlinSub2 dataset
function [ParSet] = Berlin()
    load('E:\DATA\Berlin real image\library2.mat');                         % Reference endmember bundles
    [img2d,s] = freadenvi('E:\DATA\Berlin real image\BerlinSub2');          % Original image
    load('E:\DATA\Berlin real image\RefAbundance\LAMBDA：0.000001.mat');    % Reference abundance

    ParSet.abundance = S;
    ParSet.img2d = img2d.';
    [L,ParSet.N] = size(ParSet.img2d);
    ParSet.img3d = zeros(s(1),s(2),L);
    for i = 1:L
        img3d(:,:,i) = reshape(ParSet.img2d(i,:),s(1),s(2));
    end
    ParSet.img3d = img3d;
    ParSet.endmemberBundle = library.refLevel2';
    [ParSet.row,ParSet.col,ParSet.L] = size(img3d);
    ParSet.Ptrue = 6;

    ParSet.P = 12;              % Number of clusters
    ParSet.k_means_2d = kmeans(img2d,ParSet.P,'MaxIter',400);
    ParSet.numRun = 1;          % Number of overall runs
    ParSet.LAMBDA = 0.00001;    % Parameters for colaborative sparse unmixing 
    ParSet.AL_ITERS = 2000;     % Parameters for colaborative sparse unmixing 
    ParSet.SubNum = 3;          % Number of objective functions (default setting is 3)
    ParSet.PopNum = 30;         % Population size
    ParSet.maxFE = 75000;       % Maximum number of evaluations
    ParSet.zoneL = 3;           % Parameters for Space post-processing (SSEBE)
    ParSet.percent = 0.3;       % Parameters for Space post-processing (SSEBE)
    ParSet.RMSEtf = 1;          % Whether RMSE is assessed (1 assessed, 0 not assessed)
    ParSet.th = 0;              % Spatial post-processing thresholds(Display histogram with th = 0 Manual input threshold)
    ParSet.plot = 1;            % Whether to map the evolutionary process (1 assessed, 0 not assessed)

    % Segmentation and dimensionality reduction of images based on clustering results
    % div : Set of reduced-dimensional subimages under different subsets of clusters and bands
    % divideAll ：Set of reduced-dimensional original images in different subsets of bands
    % divOrd ：Corresponding serial numbers of the pixels in the clusters to the pixels of the original image
    % len : Number of pixels in different clusters
    % imgTrans : The result of dimensionality reduction of the original image
    [ParSet.div,ParSet.divideAll,ParSet.divOrd,ParSet.len,ParSet.imgTrans] = Divide(ParSet);

    fprintf('Complete data loading and parameter setting\n');
end

%% Urban dataset
function [ParSet] = urban()
    load('E:\DATA\Urban\Urban_162.mat');
    load('E:\DATA\Urban\end6_groundTruth.mat');
    for i = 1:6
        ref{1,i} = M(:,i);
    end
    ParSet.abundance = A;
    ParSet.img2d = Y./1000;
    [L,ParSet.N] = size(ParSet.img2d);
    ParSet.img3d = zeros(nCol,nRow,nBand);
    for i = 1:L
        img3d(:,:,i) = reshape(ParSet.img2d(i,:),nCol,nRow);
    end
    ParSet.img3d = img3d;
    ParSet.endmemberBundle = ref';
    [ParSet.row,ParSet.col,ParSet.L] = size(img3d);
    ParSet.Ptrue = 6;

    ParSet.P = 12;              % Number of clusters
    ParSet.k_means_2d = kmeans(img2d,ParSet.P,'MaxIter',400);
    ParSet.numRun = 1;          % Number of overall runs
    ParSet.LAMBDA = 0.00001;    % Parameters for colaborative sparse unmixing 
    ParSet.AL_ITERS = 2000;     % Parameters for colaborative sparse unmixing 
    ParSet.SubNum = 3;          % Number of objective functions (default setting is 3)
    ParSet.PopNum = 30;         % Population size
    ParSet.maxFE = 75000;       % Maximum number of evaluations
    ParSet.zoneL = 3;           % Parameters for Space post-processing (SSEBE)
    ParSet.percent = 0.3;       % Parameters for Space post-processing (SSEBE)
    ParSet.RMSEtf = 1;          % Whether RMSE is assessed (1 assessed, 0 not assessed)
    ParSet.th = 0;              % Spatial post-processing thresholds(Display histogram with th = 0 Manual input threshold)
    ParSet.plot = 1;            % Whether to map the evolutionary process (1 assessed, 0 not assessed)
    
    %% Segmentation and dimensionality reduction of images based on clustering results
    % div : Set of reduced-dimensional subimages under different subsets of clusters and bands
    % divideAll ：Set of reduced-dimensional original images in different subsets of bands
    % divOrd ：Corresponding serial numbers of the pixels in the clusters to the pixels of the original image
    % len : Number of pixels in different clusters
    % imgTrans : The result of dimensionality reduction of the original image
    [ParSet.div,ParSet.divideAll,ParSet.divOrd,ParSet.len,ParSet.imgTrans] = Divide(ParSet);

    fprintf('Complete data loading and parameter setting\n');
end