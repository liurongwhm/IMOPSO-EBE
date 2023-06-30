function [objs] = ObjVolume(position,div)
    % Volume calculation
    [P,SubNum] = size(div);
    objs = zeros(1,SubNum);
    EndmemberReduction = zeros(6,6);
    % Calculate the value of the objective function
    for k = 1:P
        picture = div{k,1};
        EndmemberReduction(k,:) = picture(position(k),position(k+P),:);
    end
    EndmemberReduction = EndmemberReduction.';
    objs(1,1) = 1/abs(det(EndmemberReduction)/factorial(P-1));       % 1/V
end