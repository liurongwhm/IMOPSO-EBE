function [divide,divideAll,divOrd,len,imgTrans] = Divide(ParSet)
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 Pengrui Wang. You are free to use the IMOPSO-EBE for
% research purposes. All publications which use this code should reference
% "R. Liu, P. Wang, B. Du and B. Qu, "Endmember Bundle Extraction Based on 
% Improved Multi-objective Particle Swarm Optimization," in IEEE Geoscience 
% and Remote Sensing Letters, doi: 10.1109/LGRS.2023.3287919".
%--------------------------------------------------------------------------
    k_means_2d = ParSet.k_means_2d;
    img2d = ParSet.img2d;
    P = ParSet.P;
    L = ParSet.L;
    SubNum = ParSet.SubNum;
    row = ParSet.row;
    col = ParSet.col;
    N = ParSet.N;

    %% Segmentation of images in the band dimension
    subImg2d = cell(1,SubNum);
    X = cell(1,SubNum);
    [~,A] = hyperMnf(img2d, row, col);
    A = A';
    trans = A(1:P-1,:);
    SubTrans = cell(1,SubNum);
    imgTrans = [ones(1,N);trans * img2d];

    % Circular segmentation of images
    for i = 1 : SubNum
        A = (1:L);
        subImg2d{1,i} = img2d(A(i:SubNum:end),:);
        SubTrans{1,i} = trans(:,A(i:SubNum:end));
        X{1,i} = [ones(1,N);SubTrans{1,i} * subImg2d{1,i}];
    end

    %% Further separation of the segmented image by endmember category
    divide = cell(P,SubNum);
    divOrd = cell(P,1);
    for j = 1:SubNum
        for i = 1:P
            order = (1:N).';
            Xtmp = X{1,j}.';
            k_i = k_means_2d;
            k_i(k_i ~= i) = 0;
            k_i(k_i == i) = 1;
            k_tmp = k_i.*Xtmp;
            k_tmp = [k_tmp order];
            keep = k_tmp(:,1) ~= 0;
            k_tmp = k_tmp(keep,:);
            order= k_tmp(:,P+1);
            k_tmp = k_tmp(:,1:P);
            [l,~] = size(k_tmp);
            len(i) = l;
            divide{i,j} = k_tmp;
            divOrd{i,1} = order;
        end
    end

    %% Separate the whole image after dimensionality reduction by end element category
    divideAll = cell(P,1);
    for i = 1:P
        Xtmp = imgTrans.';
        k_i = k_means_2d;
        k_i(k_i ~= i) = 0;
        k_i(k_i == i) = 1;
        k_tmp = k_i.*Xtmp;
        keep = k_tmp(:,1) ~= 0;
        k_tmp = k_tmp(keep,:);
        divideAll{i,1} = k_tmp;
    end
end