function [table] = SAMpipei(A1,A2)
% Spectral angle matching :A1, A2 are the spectral matrices to be matched, A1 is the reference A2 is the estimation matrix
% size is band*p1 band*p2, p1 is the number of reference spectra, p2 is the number of estimated spectra.
[~,p1] = size(A1);
[~,p2] = size(A2);
table = zeros(4,p2+2);
for i=1:p2
    table(2,i)=i ;
    temp1=100;
    for j = 1:p1
        temp2 = SAM(A1(:,j),A2(:,i));
        if temp2<temp1
            temp1=temp2;
            table(1,i)= j;
            table(3,i)= temp2;
            table(4,i)= temp2*180/3.1415926;
        end
        
    end
end
table(3,p2+1)=sum(table(3,1:p2))/p2;
table(4,p2+1)=sum(table(4,1:p2))/p2;
table(3,p2+2)=std(table(3,1:p2));
table(4,p2+2)=std(table(4,1:p2));