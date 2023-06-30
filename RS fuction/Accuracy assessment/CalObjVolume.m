function [objs] = CalObjVolume(position,div)
    % Volume calculation
    [P,SubNum] = size(div);
    objs = zeros(1,SubNum);
    EndmemberReduction = zeros(P,P);
    % Calculate the value of the objective function
    for j = 1:SubNum
        for k = 1:P
            picture = div{k,j};
            EndmemberReduction(k,:) = picture(position(k),:);
        end
        EndmemberReduction = EndmemberReduction.';
        objs(1,j) = 1/abs(det(EndmemberReduction)/factorial(P-1));       % 1/V
    end
end