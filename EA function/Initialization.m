function [Population] = Initialization(ParSet)
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 Pengrui Wang. You are free to use the IMOPSO-EBE for
% research purposes. All publications which use this code should reference
% "R. Liu, P. Wang, B. Du and B. Qu, "Endmember Bundle Extraction Based on 
% Improved Multi-objective Particle Swarm Optimization," in IEEE Geoscience 
% and Remote Sensing Letters, doi: 10.1109/LGRS.2023.3287919".
%--------------------------------------------------------------------------
    %% Incoming parameters
    PopNum = ParSet.PopNum; 
    P = ParSet.P;
    div = ParSet.div;
    len = ParSet.len;

    %% Initializing the population
    Lower = ones(1,P);
    Upper = len;
    % Initializing empty populations
    empty_particle.position = [];                   % Particle position vector
    empty_particle.velocity = [];                   % Particle velocity vector
    empty_particle.objs = [];                       % Value of particle objective function
    empty_particle.Front = [];
    Population = repmat(empty_particle,PopNum,1);
    % Random initialization of populations
    for i = 1:PopNum
         Population(i).position = round(unifrnd(Lower,Upper));
         Population(i).velocity = zeros(1,P);
         Population(i).objs = CalObjVolume(Population(i).position,div); % Calculate the value of the objective function
         Population(i).Front = 0;
    end
end