function [indicies] = CCEBE(ParSet)
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 Pengrui Wang. You are free to use the IMOPSO-EBE for
% research purposes. All publications which use this code should reference
% "R. Liu, P. Wang, B. Du and B. Qu, "Endmember Bundle Extraction Based on 
% Improved Multi-objective Particle Swarm Optimization," in IEEE Geoscience 
% and Remote Sensing Letters, doi: 10.1109/LGRS.2023.3287919".
%--------------------------------------------------------------------------
    P = ParSet.P;
    PopNum = ParSet.PopNum;
    evaluation = ParSet.maxFE;
    evaluated = 0;
    [Population] = Initialization(ParSet);  % Initializing the population
    divOrd = ParSet.divOrd;
    IndexEndmember = zeros(PopNum,P);
    it = 1;
    %% Optimisation
    while(evaluated<evaluation)
        Offspring  = Operator(Population,ParSet);
        Population = EnvironmentalSelection([Population;Offspring],PopNum).';
        obj = cat(1,Population.objs);
    %% plot
        if ParSet.plot
            ws = ceil(sqrt(ParSet.SubNum));
            for i = 1:ParSet.SubNum
                min_list{1,i} = min(obj(:,i));
                subplot(ws,ws,i);
                record{1,i}(it) = min_list{1,i};
                plot(record{1,i});
                title(i);xlabel('Number of iterations');ylabel('<');
            end
            it = it+1;
        end
        evaluated = evaluated + size(obj,1);
    end
    decs = cat(1,Population.position);

    %% Generate 2d coordinates in the whole image according to divOrd for post-processing
    for i = 1:PopNum
        for j = 1:P
            order = divOrd{j,1};
            IndexEndmember(i,j) = order(decs(i,j)) ;
        end
    end
    indicies = unique(IndexEndmember(:));
end