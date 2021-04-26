% addpath('C:\Users\Maury-D\Documents\MATLAB\COOLPROP')
clear;
clc;
close all

div_P = 100;
% x = 0:0.01:0.1;
% x=[0, 0.1];
x = 0;
div_Q = length(x);

% FLUIDOS AVALIADOS NA FASE INICIAL
% Fluidos com baixa temperatura crítica
% fluid_l = {'R134a','R245fa','R152a','R236fa','R227ea','R143a'};
% Fluidos com média temperatura crítica
% fluid_m = {'R123','R245ca','butane','n-Pentane'};
% Fluidos com alta temperatura crítica
% fluid_h = {'Benzene','Toluene','MDM','Cyclohexane'};

% Fluidos com melhor desempenho (em termos de densidade de energia)
% fluid = {'R143a','R152a','R134a','butane','n-Pentane','Cyclohexane','Toluene','Benzene'};

fluid = {'R152a','R134a','R142b','R365mfc','R236ea','R141b'}; % TESE
fluid_std = {'R-152a','R-134a','R-142b','R-365mfc','R-236ea','R-141b'};
% fluid = {'R152a','R134a','R142b','R365mfc','R141b'}; % ARTIGO
% fluid_std = {'R-152a','R-134a','R-142b','R-365mfc','R-141b'}; % ARTIGO
% fluid_std = {'R-152a','R-134a','R-142b','R-365mfc','R236ea','R-141b'}; % TESE
% fluid = {'R152a','R134a','R142b','R365mfc','R236ea','R141b','CO2'};
% fluid = {'R134a'};
P = zeros(length(fluid),div_P);
T = zeros(length(fluid),div_P);
u = zeros(length(fluid),div_P);
D = zeros(length(fluid),div_P);
v = zeros(length(fluid),div_P);
s = zeros(length(fluid),div_P);

P_o = 101325;
T_o = 298.15;
% 
P_min = 600000;

%P_max = 3000000;
P_crit = zeros(length(fluid),1);
P_max = zeros(length(fluid),1);
u_o = zeros(length(fluid),1);
v_o = zeros(length(fluid),1);
s_o = zeros(length(fluid),1);

for k=1:length(fluid)
   P_crit(k) = CoolProp.Props1SI('Pcrit',fluid{k});
   P_max(k) = 0.9*P_crit(k);
   P(k,:) = P_min:(P_max(k)-P_min)/(div_P-1):P_max(k);
end

% ANÁLISE SUBCRÍTICA
ex = zeros(length(fluid),size(P,2),length(x)); % Exergia específica [J/kg]
Ex_kJ_kg = zeros(length(fluid),size(P,2),length(x));              % Conversão para kJ/kg 
Ex_kJ_l = zeros(length(fluid),size(P,2),length(x));
Ex_Wh_l = zeros(length(fluid),size(P,2),length(x));
Ex_Wh_kg = zeros(length(fluid),size(P,2),length(x));

for k = 1:length(fluid)
    u_o(k) = CoolProp.PropsSI('U','P',P_o,'T',T_o,fluid{k}); % Estado morto
    v_o(k) = 1/CoolProp.PropsSI('D','P',P_o,'T',T_o,fluid{k});
    s_o(k) = CoolProp.PropsSI('S','P',P_o,'T',T_o,fluid{k});

    for j=1:length(P) % Variação da pressão considerando estado de líquido saturado
        for i = 1:length(x)
            T(k,j,i) = CoolProp.PropsSI('T','P',P(k,j),'Q',x(i),fluid{k})-1;   % [K]
            u(k,j,i) = CoolProp.PropsSI('U','P',P(k,j),'Q',x(i),fluid{k});     % [J/kg]
            D(k,j,i) = CoolProp.PropsSI('D','P',P(k,j),'Q',x(i),fluid{k});     % [kg/m3]
            v(k,j,i) = 1/D(k,j,i);                                               % [m3/kg]
            s(k,j,i) = CoolProp.PropsSI('S','P',P(k,j),'Q',x(i),fluid{k});     % [J/kg-K]

            ex(k,j,i) = u(k,j,i) - u_o(k) + P_o*(v(k,j,i)-v_o(k)) - T_o*(s(k,j,i) - s_o(k)); % Exergia específica [J/kg]
            Ex_kJ_kg(k,j,i) = ex(k,j,i)./1000;              % Conversão para kJ/kg 
            Ex_kJ_l(k,j,i) = (ex(k,j,i)./1000)/(v(k,j)*1000);                % [kJ/L]
            Ex_Wh_kg(k,j,i) = ex(k,j,i)/3600;               % Conversão para Wh/kg 
            Ex_Wh_l(k,j,i) = (ex(k,j,i)/3600)/(v(k,j)*1000);                 % [Wh/L]
        end
    end
end

% CAES 
u_o_air = CoolProp.PropsSI('U','P',P_o,'T',T_o,'air'); % Estado morto
v_o_air = 1/CoolProp.PropsSI('D','P',P_o,'T',T_o,'air');
s_o_air = CoolProp.PropsSI('S','P',P_o,'T',T_o,'air');
P_CAES = P_o:125000:8000000;

T_CAES = zeros(length(P_CAES),1);    % [K]
u_CAES = zeros(length(P_CAES),1);    % [J/kg]
D_CAES = zeros(length(P_CAES),1);    % [kg/m3]
v_CAES  = zeros(length(P_CAES),1);   % [m3/kg]
s_CAES = zeros(length(P_CAES),1);    % [J/kg-K]

ex_CAES = zeros(length(P_CAES),1); % Exergia específica [J/kg]
Ex_kJ_kg_CAES = zeros(length(P_CAES),1);                                  % Conversão para kJ/kg 
Ex_kJ_l_CAES = zeros(length(P_CAES),1);                % [kJ/L]
Ex_Wh_kg_CAES = zeros(length(P_CAES),1);                                 % Conversão para Wh/kg 
Ex_Wh_l_CAES = zeros(length(P_CAES),1);                % [Wh/L]

for j=1:length(P_CAES) % Variação da pressão considerando estado de líquido saturado
    T_CAES(j) = T_o;   % [K]
    u_CAES(j) = CoolProp.PropsSI('U','P',P_CAES(j),'T',T_CAES(j),'air');     % [J/kg]
    D_CAES(j) = CoolProp.PropsSI('D','P',P_CAES(j),'T',T_CAES(j),'air');     % [kg/m3]
    v_CAES(j) = 1/D_CAES(j);                                                    % [m3/kg]
    s_CAES(j) = CoolProp.PropsSI('S','P',P_CAES(j),'T',T_CAES(j),'air');     % [J/kg-K]

    ex_CAES(j) = u_CAES(j) - u_o_air + P_o*(v_CAES(j)-v_o_air) - T_o*(s_CAES(j) - s_o_air); % Exergia específica [J/kg]
    Ex_kJ_kg_CAES(j) = ex_CAES(j)./1000;                                  % Conversão para kJ/kg 
    Ex_kJ_l_CAES(j) = (ex_CAES(j)./1000)/(v_CAES(j)*1000);                % [kJ/L]
    Ex_Wh_kg_CAES(j) = ex_CAES(j)./3600;                                   % Conversão para Wh/kg 
    Ex_Wh_l_CAES(j) = (ex_CAES(j)./3600)/(v_CAES(j)*1000);                 % [Wh/L]
end

P_LAES = P_o;
T_LAES = 77;   % [K]
u_LAES = CoolProp.PropsSI('U','P',P_LAES,'T',T_LAES,'air');     % [J/kg]
D_LAES = CoolProp.PropsSI('D','P',P_LAES,'T',T_LAES,'air');     % [kg/m3]
v_LAES = 1/D_LAES;                                              % [m3/kg]
s_LAES = CoolProp.PropsSI('S','P',P_LAES,'T',T_LAES,'air');     % [J/kg-K]

ex_LAES = u_LAES - u_o_air + P_o*(v_LAES-v_o_air) - T_o*(s_LAES - s_o_air); % Exergia específica [J/kg]
Ex_kJ_kg_LAES = ex_LAES./1000;                                  % Conversão para kJ/kg 
Ex_kJ_l_LAES = (ex_LAES./1000)/(v_LAES*1000);                % [kJ/L]
Ex_Wh_kg_LAES = ex_LAES/3600;                                   % Conversão para Wh/kg 
Ex_Wh_l_LAES = (ex_LAES/3600)/(v_LAES*1000);                 % [Wh/L]

P_kPa = P./1000;
% P_bar = P./100000;
%% Figuras gerais
for k=1:length(fluid)    % Figura da variação da densidade exergética [Wh/L] com a pressão e título
    figure('color',[1 1 1]);
    hold all
    for i=1:length(x)
        plot(P_kPa(k,:),Ex_Wh_l(k,:,i))
    end
    grid on;
    title(fluid{k})
    xlabel('P [kPa]')
    ylabel('Exergy density [Wh~L^{-1}]')
    legend(num2str(x'));
end

figure('color',[1 1 1]);
hold all
for k=1:length(fluid)   % Figura da variação da densidade exergética [Wh/L] com a pressão
    if k<=length(fluid)/2
        plot(P_kPa(k,:),Ex_Wh_l(k,:,1),'Color',k/(length(fluid)/2+1)*[1 1 1])
    else
%         plot(P_kPa(k,:),Ex_Wh_l(k,:,1),'--','Color',k/(length(fluid)+1)*[1 1 1])
        plot(P_kPa(k,:),Ex_Wh_l(k,:,1),'Color',k/(length(fluid)+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',k/(length(fluid)+1)*[1 1 1],'MarkerIndices',1:10:size(P_kPa,2))
    end
end
% plot(P_CAES./1000,Ex_Wh_l_CAES,'--k')
% plot(P_CAES./1000,ones(length(P_CAES),1)*Ex_Wh_l_LAES,'-.k')
grid on;
xtickformat('%,1.0f')
xlabel('P [kPa]')
ylabel('Exergy density [Wh~L^{-1}]')
% legenda1 = cat(2,fluid,{'Compressed Air','Liquid Air'});
% legend(legenda1)
legend(fluid_std,'Location','SouthEast')
applystyle2plot()

figure('color',[1 1 1]);
hold all
for k=1:length(fluid)   % Figura da variação da densidade exergética [Wh/kg] com a pressão
    if k<=length(fluid)/2
        plot(P_kPa(k,:),Ex_Wh_kg(k,:,1),'Color',k/(length(fluid)/2+1)*[1 1 1])
    else
%         plot(P_kPa(k,:),Ex_Wh_kg(k,:,1),'--','Color',k/(length(fluid)+1)*[1 1 1])
        plot(P_kPa(k,:),Ex_Wh_kg(k,:,1),'Color',k/(length(fluid)+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',k/(length(fluid)+1)*[1 1 1],'MarkerIndices',1:10:size(P_kPa,2))
    end
end
% plot(P_CAES./1000,Ex_Wh_kg_CAES,'--k')
% plot(P_CAES./1000,ones(length(P_CAES),1)*Ex_Wh_kg_LAES,'-.k')
grid on;
xtickformat('%,1.0f')
xlabel('P [kPa]')
ylabel('Exergy density [Wh~kg^{-1}]')
% legend(legenda1)
legend(fluid_std,'Location','SouthEast')
applystyle2plot()

figure('color',[1 1 1]);
hold all
for k=1:length(fluid)   % Figura da variação da Temperatura de saturação com a pressão
    if k<length(fluid)/2
        plot(P_kPa(k,:),T(k,:),'Color',k/(length(fluid)/2+1)*[1 1 1])
    else
        plot(P_kPa(k,:),T(k,:),'--','Color',k/(length(fluid)+1)*[1 1 1])
    end
    %plot(P_bar(i,:),T(i,:),'Color',i/(length(fluid)+1)*[1 1 1])
end
grid on;
xlabel('P [kPa]')
ylabel('T [K]') % Saturation temperature
legend(fluid)

%% Figuras gerais (kWh/m3 e kWh/kg)

figure('color',[1 1 1]);
hold all
for k=1:length(fluid)   % Figura da variação da densidade exergética [Wh/m3] com a pressão
    if k<=length(fluid)/2
        plot(P_kPa(k,:),Ex_Wh_l(k,:,1),'Color',k/(length(fluid)/2+1)*[1 1 1])
    else
%         plot(P_kPa(k,:),Ex_Wh_l(k,:,1),'--','Color',k/(length(fluid)+1)*[1 1 1])
        plot(P_kPa(k,:),Ex_Wh_l(k,:,1),'Color',k/(length(fluid)+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',k/(length(fluid)+1)*[1 1 1],'MarkerIndices',1:10:size(P_kPa,2))
    end
end
% plot(P_CAES./1000,Ex_Wh_l_CAES,'--k')
% plot(P_CAES./1000,ones(length(P_CAES),1)*Ex_Wh_l_LAES,'-.k')
grid on;
xtickformat('%,1.0f')
xlabel('P [kPa]')
ylabel('Exergy density [kWh m^{-3}]')
% legenda1 = cat(2,fluid,{'Compressed Air','Liquid Air'});
% legend(legenda1)
legend(fluid_std,'Location','SouthEast')
applystyle2plot()

figure('color',[1 1 1]);
hold all
for k=1:length(fluid)   % Figura da variação da densidade exergética [kWh/kg] com a pressão
    if k<=length(fluid)/2
        plot(P_kPa(k,:),Ex_Wh_kg(k,:,1)./1000,'Color',k/(length(fluid)/2+1)*[1 1 1])
    else
%         plot(P_kPa(k,:),Ex_Wh_kg(k,:,1),'--','Color',k/(length(fluid)+1)*[1 1 1])
        plot(P_kPa(k,:),Ex_Wh_kg(k,:,1)./1000,'Color',k/(length(fluid)+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',k/(length(fluid)+1)*[1 1 1],'MarkerIndices',1:10:size(P_kPa,2))
    end
end
% plot(P_CAES./1000,Ex_Wh_kg_CAES,'--k')
% plot(P_CAES./1000,ones(length(P_CAES),1)*Ex_Wh_kg_LAES,'-.k')
grid on;
xtickformat('%,1.0f')
xlabel('P [kPa]')
ylabel('Exergy density [kWh kg^{-1}]')
ytickformat('%,0.3f')
% legend(legenda1)
legend(fluid_std,'Location','SouthEast')
applystyle2plot()

%% Análise de potencial de custos por Wh
Fluid_Cost_per_kg = [4.49 7.01 4.41 2.00 5.00]; % Cost based on Alibaba search on 16/01/2020 - Matheus
Storage_Cost_per_m3 = 750; % Pressure tank cost [$/m3]

Ex_kWh_m3 = Ex_Wh_l;
Ex_kWh_m3_CAES = Ex_Wh_l_CAES;
Ex_kWh_kg = Ex_Wh_kg./1000;

figure('color',[1 1 1])
% subplot(1,2,1)
% subplot(2,1,1)
hold on
grid on
for k=1:length(fluid)
    Fluid_Storage_Cost_per_kWh(k,:) = Storage_Cost_per_m3./(Ex_kWh_m3(k,:,1));
    Fluid_Cost_per_kWh(k,:) = Fluid_Cost_per_kg(k)./Ex_kWh_kg(k,:,1);
    Cost(k,:) = Fluid_Storage_Cost_per_kWh(k,:) + Fluid_Cost_per_kWh(k,:);
    if k<=length(fluid)/2
        plot(P_kPa(k,:),Cost(k,:),'Color',k/(length(fluid)/2+1)*[1 1 1])
    else
        plot(P_kPa(k,:),Cost(k,:),'Color',k/(length(fluid)+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',k/(length(fluid)+1)*[1 1 1],'MarkerIndices',1:10:size(P_kPa,2))
    end
end
plot(P_CAES(:)./1000,Storage_Cost_per_m3./Ex_kWh_m3_CAES(:),'--k')
legenda = cat(2,fluid_std,{'Compressed air'});
legend(legenda)
xlabel('P [kPa]')
xtickformat('%,1.0f')
ylabel('Cost of pressure tank and fluid [USD kWh^{-1}]')
ytickformat('%,1.0f')
xlim([600 4000])
applystyle2plot()

%%
figure('color',[1 1 1])
% subplot(1,2,1)
subplot(2,1,1)
hold on
grid on
for k=1:length(fluid)
    Fluid_Storage_Cost_per_kWh(k,:) = Storage_Cost_per_m3./(Ex_kWh_m3(k,:,1));
%     Fluid_Cost_per_m3(k,:) = Fluid_Cost_per_kg(k)./v(k,:,1);
    Fluid_Cost_per_kWh(k,:) = Fluid_Cost_per_kg(k)./Ex_Wh_kg(k,:,1); 
    Cost(k,:) = Fluid_Storage_Cost_per_kWh(k,:) + Fluid_Cost_per_kWh(k,:);
    if k<=length(fluid)/2
%         plot(P_kPa(k,:),Fluid_Storage_Cost_per_m3(k,:),'Color',k/(length(fluid)/2+1)*[1 1 1])
        plot(P_kPa(k,:),Cost(k,:),'Color',k/(length(fluid)/2+1)*[1 1 1])
    else
%         plot(P_kPa(k,:),Ex_Wh_kg(k,:,1),'--','Color',k/(length(fluid)+1)*[1 1 1])
%         plot(P_kPa(k,:),Fluid_Storage_Cost_per_m3(k,:),'Color',k/(length(fluid)+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',k/(length(fluid)+1)*[1 1 1],'MarkerIndices',1:10:size(P_kPa,2))
        plot(P_kPa(k,:),Cost(k,:),'Color',k/(length(fluid)+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',k/(length(fluid)+1)*[1 1 1],'MarkerIndices',1:10:size(P_kPa,2))
    end
end
plot(P_CAES(:)./1000,Storage_Cost_per_m3./Ex_kWh_m3_CAES(:),'--k')
legenda = cat(2,fluid_std,{'Compressed air'});
legend(legenda)
xlabel('P [kPa]')
xtickformat('%,1.0f')
ylabel('Pressure tank cost [USD Wh^{-1}]')
ytickformat('%,1.1f')
applystyle2plot()

% figure('color',[1 1 1])
% subplot(1,2,2)
subplot(2,1,2)
hold on
grid on
for k=1:length(fluid)
    if k<=length(fluid)/2
        plot(P_kPa(k,:),Fluid_Cost_per_kWh(k,:),'Color',k/(length(fluid)/2+1)*[1 1 1])
    else
        plot(P_kPa(k,:),Fluid_Cost_per_kWh(k,:),'Color',k/(length(fluid)+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',k/(length(fluid)+1)*[1 1 1],'MarkerIndices',1:10:size(P_kPa,2))
    end
end
legend(fluid_std)
xlabel('P [kPa]')
xtickformat('%,1.0f')
ylabel('Fluid cost [USD~m^{-3}]')
ytickformat('%,1.1f')
applystyle2plot()

%% Análise dos componentes da exergia
i = 1;
for k=1:length(fluid)
    
    for j=1:length(P)
        ex_wh_L_u(k,j,i) = ((+ u(k,j,i) - u_o(k))/3600)/(v(k,j)*1000);
        ex_wh_L_P(k,j,i) = ((+ P_o*(v(k,j,i) - v_o(k)))/3600)/(v(k,j)*1000);
        ex_wh_L_T(k,j,i) = ((- T_o*(s(k,j,i) - s_o(k)))/3600)/(v(k,j)*1000);
    end

    figure('color',[1 1 1]) % j/kg
    plot(P(k,:,i)./1000,Ex_Wh_l(k,:,i),'k')
%     yyaxis right
    hold on
    plot(P(k,:,i)./1000,ex_wh_L_u(k,:,i),'b')
    plot(P(k,:,i)./1000,ex_wh_L_P(k,:,i),'g')
    plot(P(k,:,i)./1000,ex_wh_L_T(k,:,i),'r')
    legend('Total','u','Pv','Ts')
    grid minor
    title(fluid{k})
    xlim([800 4500])
    ylim([-100 150])
end

figure('color',[1 1 1])
hold all
for k=1:length(fluid)   % Figura da variação da densidade exergética [Wh/L] com a pressão
    if k<length(fluid)/2
        plot(P_kPa(k,:),Ex_Wh_l(k,:,1),'Color',k/(length(fluid)/2+1)*[1 1 1])
    else
        plot(P_kPa(k,:),Ex_Wh_l(k,:,1),'--','Color',k/(length(fluid)+1)*[1 1 1])
    end
end
grid minor;
xlabel('P [kPa]')
ylabel('Exergy density [Wh L^{-1}]')
legend(fluid)
%%
figure('color',[1 1 1]) % j/kg
plot((u(1,:)-u_o)./v(1,:),'b')
hold on
plot(T_o*(s(1,:)-s_o)./v(1,:),'k')
plot(P_o*(v(1,:)-v_o)./v(1,:),'r')
yyaxis right
plot((u(1,:)-u_o - T_o*(s(1,:)-s_o))./v(1,:),'g')
plot(( u(1,:)-u_o +P_o*(v(1,:)-v_o)- T_o*(s(1,:)-s_o) )./v(1,:),'k--')