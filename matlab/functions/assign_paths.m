% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function paths = assign_paths(par)
% Creates path object based on inputs:
% par.K: total number of paths
% par.type: '1d or '2d'
% 
% outputs path object with parameters
% paths.sources : total number of sources
% paths.multi : vector showing how many paths per source
% paths.AoA : angles of arrival, random on the uniform sphere
% paths.signal_vector(k).path(p) : signal vector for the p'th arrival of the k'th signal
K = par.K;

if par.forceMulti && isempty(par.forcePath)
    if K==1
        options=[1];
        r=1;
    elseif K==2
        options=[1,1];
        r=randi(1);
    elseif K==3
        options=[1,1,2;1,1,1];
        r=randi(2);
    elseif K==4
        options=[1,1,2,3;1,1,2,2;1,1,1,2;1,1,1,1];
        r=randi(4);
    elseif K==5
        options=[1,1,2,3,4;1,1,1,2,3;1,1,1,1,2;1,1,1,1,1;1,1,2,2,3;1,1,1,2,2];
        r=randi(6);
    end
elseif isempty(par.forcePath)
    if K==1
        options=[1];
        r=1;
    elseif K==2
        options=[1,2;1,1];
        r=randi(K);
    elseif K==3
        options=[1,2,3;1,1,2;1,1,1];
        r=randi(K);
    elseif K==4
        options=[1,2,3,4;1,1,2,3;1,1,2,2;1,1,1,2;1,1,1,1];
        r=randi(5);
    elseif K==5
        options=[1,2,3,4,5;1,1,2,3,4;1,1,1,2,3;1,1,1,1,2;1,1,1,1,1;1,1,2,2,3;1,1,1,2,2];
        r=randi(7);
    end
else
    options = par.forcePath;
    r=1;
end
Q=options(r,:);

paths.sources = max(Q);

for s = 1:paths.sources
    paths.multi(s) = sum(Q(:)==s);
end

% Assign azimuth and elevation (or if type = '1d', just azimuth)
switch par.type
    case '2d' 
        paths.AoA = [];
        paths.Pol = [];
        [Azi, Ele] = assignAoA(par);
        minSepFlag=0;
        
        if par.polType == "rnd_lin"
            Pol1 = pi/2 * rand(par.K, 1);
            Pol2 = zeros(par.K, 1);
        elseif par.polType == "rnd"
            Pol1 = pi/2 * rand(par.K, 1);
            Pol2 = 2*pi*rand(par.K,1)-pi;
        end
        
%         while min(diff(sort(Pol1)))<par.minSep
%             Pol1 = pi/2 * rand(par.K, 1);
%         end
%         while ~all(Azi >= par.aziRange(1) & Azi <= par.aziRange(2)) || ~all(Ele >= par.eleRange(1) & Ele <= par.eleRange(2))
%             [Azi, Ele] = assignAoA(par);
%         end
        
        if length(Azi) ~= 1
            if min(diff(sort(Azi)))<par.minSep || min(diff(sort(Ele)))<par.minSep% check if minSep criteria is met
                minSepFlag=1;
            end
        end
        
        while (~isempty(par.minSep) && minSepFlag==1) || ~all(Azi >= par.aziRange(1) & Azi <= par.aziRange(2)) || ~all(Ele >= par.eleRange(1) & Ele <= par.eleRange(2))
            [Azi, Ele] = assignAoA(par);
            if min(diff(sort(Azi)))>=par.minSep && min(diff(sort(Ele)))>=par.minSep
                minSepFlag=0;
            end
        end
        for k=1:paths.sources
            for p = 1:paths.multi(k)
                azi = Azi(1);
                ele = Ele(1);
                pol1 = Pol1(1);
                pol2 = Pol2(1);
                paths.AoA = [paths.AoA; [ele, azi]];
                paths.Pol = [paths.Pol; [pol1, pol2]];
                [~,~,paths.signal_vector(k).path(p,:)] = VectorSensor([ele,azi], [pol1, pol2]);
                Azi(1)=[];
                Ele(1)=[];
                Pol1(1) = [];
                Pol2(1) = [];
            end               
        end     
    case '1d'
        paths.AoA = [];
        paths.Pol = [];
        Azi = assignAzi(par);
        minSepFlag=0;
        
        if par.polType == "rnd_lin"
            Pol1 = pi/2 * rand(par.K, 1);
            Pol2 = zeros(par.K, 1);
        elseif par.polType == "rnd"
            Pol1 = pi/2 * rand(par.K, 1);
            Pol2 = 2*pi*rand(par.K,1)-pi;
        end
        
%         while 
%             Azi = assignAzi(par);
%         end
        
        if min(diff(sort(Azi)))<par.minSep % check if minSep criteria is met
            minSepFlag=1;
        end
        while (~isempty(par.minSep) && minSepFlag==1) || ~all(Azi >= par.aziRange(1) & Azi <= par.aziRange(2))
            Azi = assignAzi(par);
            if min(diff(sort(Azi)))>=par.minSep
                minSepFlag=0;
            end
        end
        for k=1:paths.sources
            for p = 1:paths.multi(k)
                azi = Azi(1);
                ele = pi/2;
                pol1 = Pol1(1);
                pol2 = Pol2(1);
                paths.AoA = [paths.AoA; [ele, azi]];
                paths.Pol = [paths.Pol; [pol1, pol2]];
                [~,~,paths.signal_vector(k).path(p,:)] = VectorSensor([ele,azi], [pol1, pol2]);
                Azi(1)=[];
                Pol1(1) = [];
                Pol2(1) = [];
            end               
        end        
end



end

function azi = assignAzi(par)
    if isempty(par.forcePath)
        k = par.K;
    else
        k = length(par.forcePath);
    end
    for i=1:k
        azi(i) = rand*2*pi;
    end
end

function [azi, ele] = assignAoA(par)
    if isempty(par.forcePath)
        k = par.K;
    else
        k = length(par.forcePath);
    end
    vector = randn(3, k);
    unitVector = vector./sqrt(sum(vector.^2,1));
    [azi,ele] = cart2sph(unitVector(1,:),unitVector(2,:),unitVector(3,:));
    ele = ele + pi/2;    
    azi = azi + pi;
end