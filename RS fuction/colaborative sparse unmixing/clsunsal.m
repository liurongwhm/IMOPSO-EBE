function [x,res_p,res_d] = clsunsal(M,y,varargin)

%% [x] = clsunsal_v1(M,y,varargin)
%
%  CLSUNSAL -> collaborative sparse unmixing via variable splitting and augmented
%  Lagrangian  
%
%% --------------- Description --------------------------------------------
%
%  CLSUNSAL solves the following l2-l1 optimization  problem 
%  [size(M) = (L,p); size(X) = (p,N)]; size(Y) = (L,N)]
%
%         min  (1/2) ||M X-Y||^2_F + lambda ||X||_{2,1}
%          X              
%
%  where ||X||_{2,1} = sum(sqrt(sum(X.^2),2))
% 
%    CONSTRAINTS ACCEPTED:
%
%    1) POSITIVITY:  X >= 0;
%    2) ADDONE:  sum(X) = ones(1,N);
%
%          
%
%% -------------------- Line of Attack  -----------------------------------
%
%  CLSUNSAL solves the above optimization problem by introducing a variable
%  splitting and then solving the resulting constrained optimization with
%  the augmented Lagrangian method of multipliers (ADMM). 
% 
% 
%         min  (1/2) ||M X-Y||^2_F + lambda ||V1||_{2,1} + i_+(V2)
%          X,Z              
%         subject to: sum(X) = ones(1,N)); V1 = X,  V2 = X;
%
%  Augmented Lagrangian (scaled version):
%
%       L(X,V1,V2,D1,D2) = (1/2) ||M X-Y||^2_F 
%                      + lambda ||V1||_{2,1} + mu/2||X-V1-D1||^2_F
%                      + i_+(V2) + mu/2||X-V2-D2||^2_F
%       
%  where D1 and D2 are the scaled Lagrange multipliers
%
%
%  ADMM:
%
%      do 
%        X  <-- arg min L(X,V1,V2,D1,D2)
%                    X,   s.t: sum(X) = ones(1,N));
%        V1  <-- arg min L(X,V1,V2,D1,D2)
%                     V1
%        V2  <-- arg min L(X,V1,V2,D1,D2)
%                     V2
%        D1  <-- D1 - (X-V1);
%        D2  <-- D2 - (X-V2);
%
%      while ~stop_rulde
%  
% For details see
%
% M. D. Iordache, J. M. Bioucas-Dias and A. Plaza. "Collaborative Sparse 
% Regression for Hyperspectral Unmixing",  IEEE Transactions on Geoscience 
% and Remote Sensing, 2013  (accepted).
% http://dx.doi.org/10.1109/TGRS.2013.2240001 
%
% NOTE: this version  differs from that  presented in the paper in two
% points: 
%     
%    a) the sum-to-one constraint is allowed
%
%    b) use just two spliting instead of the three proposed in the paper.
%       
% The version clsunsal_v1 implemts exactly the splitting scheme shwon in
% the paper
%
% ------------------------------------------------------------------------
%%  ===== Required inputs =============
%
%  M - [L(channels) x p(endmembers)] mixing matrix
%
%  y - matrix with  L(channels) x N(pixels).
%      each pixel is a linear mixture of p endmembers
%      signatures y = M*x + noise,
%
%      
%
%
%%  ====================== Optional inputs =============================
%
%  'AL_ITERS' - Minimum number of augmented Lagrangian iterations
%               Default: 100;
%               
%  lambda - regularization parameter. 
%
%
%  'POSITIVITY'  = {'yes', 'no'}; Enforces the positivity constraint: 
%                   X >= 0
%                   Default 'no'
%
%  'ADDONE'  = {'yes', 'no'}; Enforces the positivity constraint: X >= 0
%              Default 'no'
% 
%   'TOL'    - tolerance for the primal and  dual residuals 
%              Default = 1e-4; 
%
%
%  'verbose'   = {'yes', 'no'}; 
%                 'no' - work silently
%                 'yes' - display warnings
%                  Default 'no'
%        
%%  =========================== Outputs ==================================
%
% X  =  [pxN] estimated mixing matrix
%
%

%%
% ------------------------------------------------------------------
% Author: Jose Bioucas-Dias, 2012 
%
% 
% Modification: 
%
%     Jose Bioucas-Dias, July 2014 - corection of the projection onto R_+
%
%     

%% -------------------------------------------------------------------------
%
% Copyright (July, 2012):        Jos?Bioucas-Dias (bioucas@lx.it.pt)
%
% CLSUNSAL is distributed under the terms of
% the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ---------------------------------------------------------------------



%%
%--------------------------------------------------------------
% test for number of required parametres
%--------------------------------------------------------------
if (nargin-length(varargin)) ~= 2
    error('Wrong number of required parameters');
end
% mixing matrixsize
[LM,p] = size(M);
% data set size
[L,N] = size(y);
if (LM ~= L)
    error('mixing matrix M and data set y are inconsistent');
end
% if (L<p)
%     error('Insufficient number of columns in y');
% end


%%
%--------------------------------------------------------------
% Set the defaults for the optional parameters
%--------------------------------------------------------------
% maximum number of AL iteration
AL_iters = 1000;
% regularizatio parameter
lambda = 0.0;
% display only sunsal warnings
verbose = 'off';
% Positivity constraint
positivity = 'no';
% Sum-to-one constraint
addone = 'no';
% tolerance for the primal and dual residues
tol = 1e-4;
% initialization
x0 = 0;

%%
%--------------------------------------------------------------
% Local variables
%--------------------------------------------------------------


%--------------------------------------------------------------
% Read the optional parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'AL_ITERS'
                AL_iters = round(varargin{i+1});
                if (AL_iters <= 0 )
                       error('AL_iters must a positive integer');
                end
            case 'LAMBDA'
                lambda = varargin{i+1};
                if lambda < 0
                       error('lambda must be positive');
                end
            case 'POSITIVITY'
                positivity = varargin{i+1};
            case 'ADDONE'
                addone = varargin{i+1};
            case 'TOL'
                tol = varargin{i+1};
            case 'VERBOSE'
                verbose = varargin{i+1};
            case 'X0'
                x0 = varargin{i+1};
                if (size(x0,1) ~= p) || (size(x0,1) ~= N)
                    error('initial X is  inconsistent with M or Y');
                end
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end


% compute mean norm
norm_y = sqrt(mean(mean(y.^2)));
% rescale M and Y and lambda
M = M/norm_y;
y = y/norm_y;
lambda = lambda/norm_y^2;

  

%%
%---------------------------------------------
% just least squares
%---------------------------------------------
if sum(sum(lambda == 0)) &&  strcmp(positivity,'no') && strcmp(addone,'no')
    z = pinv(M)*y;
    % primal and dual residues
    res_p = 0;
    res_d = 0;
    return
end
%---------------------------------------------
% least squares constrained (sum(x) = 1)
%---------------------------------------------
SMALL = 1e-12;
B = ones(1,p);
a = ones(1,N);

if  strcmp(addone,'yes') && strcmp(positivity,'no') 
    F = M'*M;
    % test if F is invertible
    if rcond(F) > SMALL
        % compute the solution explicitly
        IF = inv(F);
        z = IF*M'*y-IF*B'*inv(B*IF*B')*(B*IF*M'*y-a);
        % primal and dual residues
        res_p = 0;
        res_d = 0;
        return
    end
end


%%
%---------------------------------------------
%  Constants and initializations
%---------------------------------------------
mu_AL = 0.01;
mu = 10*mean(lambda(:)) + mu_AL;

%F = M'*M+mu*eye(p);
[UF,SF] = svd(M'*M);
sF = diag(SF);
IF = UF*diag(1./(sF+2*mu))*UF';
%IF = inv(F);
Aux = IF*B'*inv(B*IF*B');
x_aux = Aux*a;
IF1 = (IF-Aux*B*IF);


yy = M'*y;

%%
%---------------------------------------------
%  Initializations
%---------------------------------------------

% no intial solution supplied
if x0 == 0
    x= IF*M'*y;
end

% auxiliary variables
v1 = x;
v2 = x;

% scaled Lagrange Multipliers
d1  = 0*v1;
d2  = 0*v2;

%%
%---------------------------------------------
%  AL iterations - main body
%---------------------------------------------
tol1 = sqrt(N*p)*tol;
tol2 = sqrt(N*p)*tol;
i=1;
res_p = inf;
res_d = inf;
mu_changed = 0;
while (i <= AL_iters) && ((abs (res_p) > tol1) || (abs (res_d) > tol2)) 
    % save z to be used later
    if mod(i,10) == 1
        v10 = v1;
        v20 = v2;
    end
    % minimize with respect to v1
    v1 =  vector_soft_row(x-d1,lambda/mu);
    % minimize wrt v2
    v2 = x-d2;
    % teste for positivity
    if strcmp(positivity,'yes')
        v2 = max(0,v2);
    end
    % teste for sum-to-one 
    if strcmp(addone,'yes')
       x = IF1*(yy + mu*(v1+d1)+ mu*(v2+d2))+x_aux;
    else
       x = IF*(yy + mu*(v1+d1)+ mu*(v2+d2));
    end
    
    % Lagrange multipliers update
    d1 = d1 -(x-v1);
    d2 = d2 -(x-v2);

    % update mu so to keep primal and dual residuals whithin a factor of 10
    if mod(i,10) == 1
        % primal residue
        res_p = sqrt(norm(x-v1,'fro')^2 + norm(x-v2,'fro')^2);
        % dual residue
        res_d = mu*norm(v1-v10+v2-v20,'fro');
        if  strcmp(verbose,'yes')
            fprintf(' i = %f, res_p = %f, res_d = %f\n',i,res_p,res_d)
        end
        % update mu
        if res_p > 10*res_d
            mu = mu*2;
            d1 = d1/2;
            d2 = d2/2;
            mu_changed = 1;
        elseif res_d > 10*res_p
            mu = mu/2;
            d1 = d1*2;
            d2 = d2*2;
            mu_changed = 1;
        end
        if  mu_changed
           % update IF and IF1

           IF = UF*diag(1./(sF+2*mu))*UF';
           Aux = IF*B'*inv(B*IF*B');
           x_aux = Aux*a;
           IF1 = (IF-Aux*B*IF);
            mu_changed = 0;
            
        end
        
        
    end
    
    i=i+1;
        
   
       
end

    
 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
