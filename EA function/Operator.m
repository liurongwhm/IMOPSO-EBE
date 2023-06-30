function Offspring = Operator(Population,ParSet,Parameter)
% The particle swarm optimization in IMOPSO-EBE
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 Pengrui Wang. You are free to use the IMOPSO-EBE for
% research purposes. All publications which use this code should reference
% "R. Liu, P. Wang, B. Du and B. Qu, "Endmember Bundle Extraction Based on 
% Improved Multi-objective Particle Swarm Optimization," in IEEE Geoscience 
% and Remote Sensing Letters, doi: 10.1109/LGRS.2023.3287919".
%--------------------------------------------------------------------------
    if nargin > 2
        [proM,disM] = deal(Parameter{:});
    else
        [proM,disM] = deal(1,20);
    end
    objs = cat(1,Population.objs);
    D = ParSet.SubNum;            
    P = ParSet.P;
    PopNum = ParSet.PopNum;
    len = ParSet.len;
    OffPopNum = 2*PopNum;
    PlNobjs = zeros(1,2);   Off_V = zeros(OffPopNum,P);   Off_P = zeros(OffPopNum,P);
    %% Get leaders  
    Front = NDSort(objs,inf);
    for i = 1:PopNum
        Population(i).Front = Front(1,i);
    end    
    
    %% Competition
    for i = 1 : PopNum
        r1 = rand(1,P);
        r2 = rand(1,P);
        PlN = randperm(PopNum,2);
        if Population(PlN(1)).Front ~= Population(PlN(2)).Front
            mask   = (Population(PlN(1)).Front < Population(PlN(2)).Front);
            winner = mask.*PlN(1) + ~mask.*PlN(2);  loser  = ~mask.*PlN(1) + mask.*PlN(2);
        else
            for j = 1:2
            PlNobjs(j) = CalObjVolume(Population(PlN(j)).position,ParSet.divideAll);
            end
            mask   = (PlNobjs(1) > PlNobjs(2));
            winner = ~mask.*PlN(1) + mask.*PlN(2);  loser  = ~mask.*PlN(2) + mask.*PlN(1);
        end
        Off_V(i,:) = r1.*Population(loser).velocity + r2.*(Population(winner).position-Population(loser).position);
        Off_P(i,:) = Population(loser).position + Off_V(i,:);
    end
    
    %% Cooperation 1
    Site  = rand(PopNum,P) < proM/D;
    mu    = rand(PopNum,P);
    temp  = Site & mu<=0.5;
    rowrank = randperm(PopNum);
    Off_P(i+1:i+PopNum,:) = Off_P(1:PopNum,:).*(1-temp) + Off_P(rowrank,:).*temp;
    %% Cooperation 2
%     for i = 1 : ceil(PopNum)
%         a = round(rand(1,P));
%         PlN = randperm(PopNum,2);
%         Off_P(PopNum+2*i,:) = Population(PlN(1)).position.*a + Population(PlN(2)).position.*(1-a);
%         Off_P(PopNum+2*i-1,:) = Population(PlN(2)).position.*a + Population(PlN(1)).position.*(1-a);
%         Off_V(PopNum+2*i,:) = Population(PlN(1)).velocity.*a + Population(PlN(2)).velocity.*(1-a);
%         Off_V(PopNum+2*i-1,:) = Population(PlN(2)).velocity.*a + Population(PlN(1)).velocity.*(1-a);
%     end

    %% Polynomial mutation
    Lower = repmat(ones(1,P),OffPopNum,1);
    Upper = repmat(len,OffPopNum,1);
    Off_P = min(max(Off_P,Lower),Upper);

    Site  = rand(OffPopNum,P) < proM/D;
    mu    = rand(OffPopNum,P);
    temp  = Site & mu<=0.5;
    Off_P(temp) = Off_P(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                  (1-(Off_P(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Off_P(temp) = Off_P(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                  (1-(Upper(temp)-Off_P(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
	Offspring = SOLUTION(Off_P,Off_V,ParSet);
end