function [indicies] = SSEBE( img3d,w,percent,idx_candidate,th)
%An Image-Based Endmember Bundle Extraction Algorithm Using Both Spatial and Spectral Information
% input: img3d
%w the size of blocks
%percent: the percentage of endmember numbers from each block
%in the original paper
% blocksize of synthetic image is 55, threshold 0.6, projection number of PPI is 10000
%cuprite data: blocksize 30, threshold 0.01, PPI 10000
%percent is 0.03

w1 = w;
w2 = w;
[m, n, k]=size(img3d);
img=zeros(1,m*n);
img(idx_candidate)=1;
[mark] = hyperConvert3d(img,m, n, 1);
% imshow(mark);

%% Calculating HI
HI_SID=ones(m,n)*500;  % initial value:Homogeneity index
d=fix(w1/2);
count=1;

for i=d+1:m-d
    for j=d+1:n-d
        if mark(i,j)==1
            index=1;
            center=reshape(img3d(i,j,:), k,1);
            for i1=-d:d
                for j1=-d:d
                    neighbor=reshape(img3d(i+i1,j+j1,:), k,1);
                    calSID(index)  = real(SID( center , neighbor ));
                    index=index+1;
                end
            end
            HI_SID(i,j)=max(calSID);
            testSID(count)=max(calSID);
            count=count+1;
        end
    end
end
testSID(testSID>500)=[];

%% Set thresholds based on HI histogram
if th ~= 0
    threshold0 =th;
else
    figure;
    title('HI threshold selection');xlabel('Number of endmembers');ylabel('HI');
    histogram(testSID,20);
    fprintf('Please enter tolerance\n');
    threshold0 = input(' ');
end
indicies=zeros(m,n);
count=1;

%% Space handling
for i=w2+d+1:w2:m+w2-1-d
    for j=w2+d+1:w2:n+w2-1-d
        if ( i<= m && j <= n )
            height_S=w2;
            width_S=w2;
        else if(i<= m && j > n)
                height_S=w2;
                width_S=w2-(j-n);
            else if ( i > m && j <= n )
                    height_S=w2-(i-m);
                    width_S=w2;
                else if (i> m && j > n )
                        height_S=w2-(i-m);
                        width_S=w2-(j-n);
                    end
                end
            end
        end
        
        threshold=threshold0;
        temp=size((find (HI_SID(i+1-w2:i+height_S-w2,j+1-w2:j+width_S-w2)<threshold)),1);
        while (temp > percent* height_S*width_S)
            threshold=0.9*threshold;
            asd = HI_SID(i+1-w2:i+height_S-w2,j+1-w2:j+width_S-w2);
            temp=size((find (asd<threshold & asd>0)),1);
        end
        for mm=1-w2:-w2+height_S
            for nn=1-w2:-w2+width_S
                if (HI_SID(i+mm,j+nn) < threshold)
                    temp2=reshape(img3d(i+mm,j+nn,:), k,1);
                    spec(:,count)=temp2;
                    indicies(i+mm,j+nn)=1;
                    count=count+1;
                end
            end
        end
    end
end
indicies = reshape(indicies,1,m*n);
indicies = find(indicies==1);
end


% hyperPpi
function [val, idx] = hyperPpi(M, numSkewers)
[p, N] = size(M);
% Remove data mean
u = mean(M.').';
M = M - repmat(u, 1, N);
% Generate skewers
skewers = randn(p, numSkewers);
votes = zeros(N, 1);
for kk=1:numSkewers
    % Project all the data onto a skewer
    tmp = abs(skewers(:,kk).'*M);
    [val, idx] = max(tmp);
    votes(idx) = votes(idx) + 1;
end
[val, idx] = sort(votes, 'descend');
end

% SID
function [ SID ] = SID( a , b )
if length(a) ~= length(b)
    error('the dimensionality of two inputs are not equal');
end
p=a./sum(a);
q=b./sum(b);
%    p=a;
%    q=b;
SID = sum((p-q).*log(p./q));
end

% hyperConvert3d
function [img] = hyperConvert3d(img, h, w, numBands)
if (ndims(img) ~= 2)
    error('Input image must be p x N.');
end
[numBands, N] = size(img);
if (1 == N)
    img = reshape(img, h, w);
else
    img = reshape(img.', h, w, numBands);
end
return;
end
