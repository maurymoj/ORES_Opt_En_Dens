clear;
clc;
close all;

% FLUIDOS AVALIADOS NA FASE INICIAL
% Fluidos com baixa temperatura crítica
% fluid_l = {'R134a','R245fa','R152a','R236fa','R227ea','R143a'};
% Fluidos com média temperatura crítica
% fluid_m = {'R123','R245ca','butane','n-Pentane'};
% Fluidos com alta temperatura crítica
% fluid_h = {'Benzene','Toluene','MDM','Cyclohexane'};

% Fluidos com melhor desempenho (em termos de densidade de energia)
% fluid = {'R143a','R152a','R134a','butane','n-Pentane','Cyclohexane','Toluene','Benzene'};

fluid = {'R152a','R134a','R142b','R365mfc','R236ea','R141b'}; % ARTIGO

% Properties
Properties = zeros(length(fluid),7);

for i=1:length(fluid)
    Properties(i,1) = 1000*CoolProp.Props1SI(fluid{i},'M'); % Molar mass
    Properties(i,2) = CoolProp.PropsSI('T','P',101325,'Q',0,fluid{i});                                   % 
    Properties(i,3) = CoolProp.Props1SI(fluid{i},'Tcrit');  % Critical T
    Properties(i,4) = CoolProp.Props1SI(fluid{i},'Pcrit');  % Critical P
    Properties(i,5) = 0;                                    %
    Properties(i,6) = 0;%CoolProp.Props1SI(fluid{i},'ODP'); % Ozone Depletion Potential
    Properties(i,7) = CoolProp.Props1SI(fluid{i},'GWP100'); % Global Warming Potential (100years) 
end