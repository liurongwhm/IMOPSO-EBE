function Offspring = SOLUTION(PopPos,AddPro,ParSet)
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 Pengrui Wang. You are free to use the IMOPSO-EBE for
% research purposes. All publications which use this code should reference
% "R. Liu, P. Wang, B. Du and B. Qu, "Endmember Bundle Extraction Based on 
% Improved Multi-objective Particle Swarm Optimization," in IEEE Geoscience 
% and Remote Sensing Letters, doi: 10.1109/LGRS.2023.3287919".
%--------------------------------------------------------------------------
    PopNum = 2*ParSet.PopNum;
    Upper = ParSet.len;
    if nargin > 0
        empty_particle.position = [];   % Particle position vector
        empty_particle.velocity = [];   % Particle velocity vector
        empty_particle.objs = [];       % Value of particle objective function
        Offspring = repmat(empty_particle,PopNum,1);
        PopPos = floor(PopPos);
        PopPos = max(min(PopPos,Upper),ones(size(PopPos,1),1));
        
        % 计算目标值
        PopObj = zeros(PopNum,ParSet.SubNum);
        for i = 1 : PopNum
            PopObj(i,:) = CalObjVolume(PopPos(i,:),ParSet.div);     % Calculate the value of the objective function
        end
        Front = NDSort(PopObj,inf);
        for i = 1 : PopNum
            Offspring(i).position   = PopPos(i,:);
            Offspring(i).objs       = PopObj(i,:);
            Offspring(i).Front      = Front(1,i);
            Offspring(i).velocity   = AddPro(i,:);
        end
    end
end