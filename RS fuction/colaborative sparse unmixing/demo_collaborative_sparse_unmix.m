%  This demo illusttates the collabotative sparse regression concepts 
%  explored in the  paper:
%
%
%  M. D. Iordache, J. M. Bioucas-Dias and A. Plaza. "Collaborative Sparse 
%  Regression for Hyperspectral Unmixing",  IEEE Transactions on Geoscience 
%  and Remote Sensing, 2013  (accepted). 
%  http://dx.doi.org/10.1109/TGRS.2013.2240001 
%
%  -----------------------------------------------------------------------
%  
%  Colaborative spase regression problem:
%
%  Let  -> Y [L,np]  matrix containing the observed spectral vectors i its
%                   columns
%
%       -> A  [L,m] mixing matrix with K spectral signatures (usually m > L)
% 
%       -> X  [m,np] abundance matrix 
%
%
%  Optimization problem: 
%
%   min  (1/2)  ||AX-Y||^2_F + lambda ||X||_{2,1}
%     X
%
%    OPTIONAL CONSTRAINTS:
%
%    1) Positivity X >= 0;
%    2) Sum-To-One sum(X) = 1
%
%   where  ||X||_{2,1} = sum(sqrt(sum(X.^2,2))) is the l_{2,1} mixed 
%   norm with promotes   rows of zeros in X.
%   
%
% Author: Jose Bioucas-Dias, February, 2013
%

clear all
close all

% define random states for  reproducible result
rand('state',23);
randn('state',23);
%


%--------------------------------------------------------------------------
% Load Library (matrix A)
%--------------------------------------------------------------------------
% 1 - USGS           (L = 224; m = 498)
% 2 - USGS - pruned  (L = 224; m = 342) (3 deg)
% 3 - USGS - pruned  (L = 38;  m = 62) (10 deg)
% 4 - USGS - pruned  (L = 8;   m = 12) (20 deg)
% 5 - USGS - pruned  (L = 224; m = 60) (3 deg)
% 6 - iid Gaussian
% 7 - iid Uniform

% size of the mixing matrix [Lxm] (only applies to libraries 6 and 7) 
% undetermined matrix
L = 200;
m = 400;

% set library
library = 3;    %
% number of pixels
np = 100;
% sparsity of the source
p = 10;

% Y = AX + N;
% set SNR in dBs  (SNR = ||MX||^2_F / ||N||^2_F)
SNR = 60; 

% 
SHAPE_PARAMETER = 1;      %(Dirichlet parameter) abundances uniformely 
                          %distribted over the simplex
MAX_PURIRY = 1;           % do not theshold abundances
OUTLIERS   = 0;           % no outliers
PURE_PIXELS = 'no';       % no pure pixels



%% --------------------------------------------------------------------------
%       Start simulation
%--------------------------------------------------------------------------

switch library
    case 1  % A_1
        load USGS_1995_Library.mat
        wavelengths = datalib(:,1);
        [dummy, indexes] = sort(wavelengths);
        A = datalib(indexes,4:end);
        names = names(4:end,:);
        clear datalib;
    case 2
        load USGS_pruned_3_deg.mat
        A = B;
        clear B;
    case 3
         load USGS_pruned_10_deg.mat
          A = B(1:6:end,:);
         clear B;
    case 4
         load USGS_pruned_20_deg.mat
         A = B(1:30:end,:);
         clear B;
    case 5
         load USGS_pruned_30_deg.mat
         A = B;
        clear B;
    case 6
        A = randn(L,2*L);
    case 7
        A = rand(L,2*L);
    otherwise
        disp('Unknown library')
end

[L,m] = size(A);  % L = number of bands; m = number of materials
% normalize A
% nA = sqrt(sum(A.^2));
% A = A./repmat(nA,L,1);


%% 
% -------------------------------------------------------------------
%            Generate data 
% ------------------------------------------------------------------- 

% mixing matrix
index = randperm(m);
M = A(:,index(1:p));
   
% generate the data
[Y,Xaux,N] = spectMixGen(M,np, ...
         'Source_pdf', 'Diri_id', ...
         'pdf_pars',SHAPE_PARAMETER,...
         'max_purity',MAX_PURIRY*ones(1,p), ...
         'no_outliers',OUTLIERS, ...
         'pure_pixels', PURE_PIXELS, ...
         'violation_extremes',[1,1.2], ....
         'snr', SNR, ...
         'noise_shape','uniform');
     
 % write X wrt A
 X = zeros(m,np);
 X(index(1:p),:) = Xaux;
         
%%
%--------------------------------------------------------------------------
%       Project data on the affine set defined by the data in the l2 sense
%-------------------------------------------------------------------------
%
%   The application of this projection ensures that the data belongs to
%   an affine set.
%
%   Up is an isometric matrix that spans the subspace where Y lives
%   Yp contains the coordinates wrt Up

%[Yp,Up,my,angles,scales] = affineProj(Y,p,'proj_type','orth');


%%
%--------------------------------------------------------------------------
% SOLVERS
%--------------------------------------------------------------------------

% least squares
[X_hat] = sunsal(A,Y);


% positivity, addone, ||X||_1 (component-wise sparsity)
[X_hat_l11] = sunsal(A,Y,'POSITIVITY','yes','VERBOSE','yes','ADDONE','no', ...
    'lambda', 1e-4,'AL_ITERS',2000, 'TOL', 1e-6);

% positivity, addone, ||X||_{2,1} (collaborative sparsity)
[X_hat_l21] = clsunsal(A,Y,'POSITIVITY','yes','VERBOSE','yes','ADDONE','yes', ...
    'lambda', 3e-4,'AL_ITERS',2000, 'TOL', 1e-8);



%%
%--------------------------------------------------------------------------
% Print Results
%--------------------------------------------------------------------------
fprintf('\nLS_ERR: \nGSUNSAL = %f\n', norm(X_hat-X,1));
fprintf('\nL12_ERR: \nGSUNSAL_{l11} = %f\n', norm(X_hat_l11-X,1));
fprintf('\nL12_ERR: \nGSUNSAL_{l22} = %f\n', norm(X_hat_l21-X,1));

xx = 1:m;
figure(1);

subplot(131)
imagesc(X)
title('Original X')

subplot(132)
imagesc(X_hat_l11)
title('Component-wise sparsity')


subplot(133)
imagesc(X_hat_l21)
title('row sparsity')

% plot sum of X_hat
figure(2);
subplot(121)
hold on
stem(1:m, sum(abs(X),2)','b')
stem(1:m, sum(abs(X_hat_l11),2)','r')
title('sum(X-hat,2)  SUNSAL_{l11}')

subplot(122)
hold on
stem(1:m, sum(abs(X),2)','b')
stem(1:m, sum(abs(X_hat_l21),2)','r')
title('sum(X-hat,2)  GLSUNSAL_{l21}')

% support detection

figure(3)
hold on
support = zeros(1,m);
support(index(1:p)) = 3;

stem(1:m, support,'b')
stem(1:m, ((sum(abs(X_hat_l21),2)' > np/p/20))*2, 'r')
stem(1:m, (sum(abs(X_hat_l11),2)' > np/p/20), 'g')
title('Support detection')
legend('X','l21', 'l1')





