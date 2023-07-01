load end4_groundTruth.mat;

nEnd = length (cood);

for i = 1 : nEnd
    subplot (2,4, i);
    plot (M(:,i));
    title(cood(i));
    grid on;
    
    a2d = reshape (A(i,:), [nRow nCol]);
    subplot (2,4, i+nEnd);
    imshow (a2d, []);
    title(cood(i));
end