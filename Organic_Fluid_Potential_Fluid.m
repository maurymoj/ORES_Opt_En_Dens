% Property table
clear;
clc;
close all

%fluid = {'R152a','R134a','R142b','R365mfc','R236ea','R141b'}; % TESE
%fluid_std = {'R-152a','R-134a','R-142b','R-365mfc','R-236ea','R-141b'};
fluid = {'R152a','R134a','R142b','R365mfc','R141b'}; % TESE
fluid_std = {'R-152a','R-134a','R-142b','R-365mfc','R-141b'};

for i=1:length(fluid)
    f = fluid{i};
    M_m = 1000*CoolProp.Props1SI(f,'M');
    T_c = CoolProp.Props1SI(f,'Tcrit');
    P_c = CoolProp.Props1SI(f,'Pcrit')/1000;
%     ODP = CoolProp.Props1SI(f,'ODP');
    GWP = CoolProp.Props1SI(f,'GWP100');
%     {f, M_m, T_c,P_c,GWP}
end

%% Exergy density
% clc
n_fluids = length(fluid);

x=0:0.1:0.5;
Volume = 2:2:10;%[1, 2, 3, 5, 8]; [m3]

div_P = 100;
div_Q = length(x);
div_Vol = length(Volume);

% Inicialização das variáveis
P_crit = zeros(n_fluids,1);
P_max = zeros(n_fluids,1);
P_min = zeros(n_fluids,1);
u_o = zeros(n_fluids,1);
v_o = zeros(n_fluids,1);
s_o = zeros(n_fluids,1);

P = zeros(n_fluids,div_P);
T = zeros(n_fluids,div_P,div_Q);
u = zeros(n_fluids,div_P,div_Q);
D = zeros(n_fluids,div_P,div_Q);
v = zeros(n_fluids,div_P,div_Q);
s = zeros(n_fluids,div_P,div_Q);

C_tank = zeros(n_fluids,div_P,div_Q,div_Vol);
m_tank = zeros(n_fluids,div_P,div_Q,div_Vol);
C_fluid = zeros(n_fluids,div_P,div_Q,div_Vol);
C_st = zeros(n_fluids,div_P,div_Q,div_Vol);
Ex_st = zeros(n_fluids,div_P,div_Q,div_Vol);

% P_CAES = zeros(n_fluids,div_P);
T_CAES = zeros(n_fluids,div_P);
u_CAES = zeros(n_fluids,div_P);
D_CAES = zeros(n_fluids,div_P);
v_CAES = zeros(n_fluids,div_P);
s_CAES = zeros(n_fluids,div_P);

ex_u = zeros(n_fluids,div_P,div_Q);
ex_v = zeros(n_fluids,div_P,div_Q);
ex_s = zeros(n_fluids,div_P,div_Q);
ex = zeros(n_fluids,div_P,div_Q); % Exergia específica [J/kg]
Ex_kJ_kg = zeros(n_fluids,div_P,div_Q);              % Conversão para kJ/kg 
Ex_kJ_l = zeros(n_fluids,div_P,div_Q);
Ex_Wh_l = zeros(n_fluids,div_P,div_Q);
Ex_Wh_kg = zeros(n_fluids,div_P,div_Q);

ex_CAES = zeros(div_P,1); % Exergia específica [J/kg]
Ex_kJ_kg_CAES = zeros(div_P,1);               % Conversão para kJ/kg 
Ex_kJ_l_CAES = zeros(div_P,1);                % [kJ/L]
Ex_Wh_kg_CAES = zeros(div_P,1);               % Conversão para Wh/kg 
Ex_Wh_l_CAES = zeros(div_P,1);                % [Wh/L]

C_st_CAES = zeros(div_P,div_Vol);
Ex_st_CAES = zeros(div_P,div_Vol);

% Estado morto
P_o = 101325;
T_o = 298.15;
%-------------------------- ANÁLISE SUBCRÍTICA ---------------------------%

% Definição dos limites de pressão e vetor de valores de pressão para cada
% fluido
for k=1:n_fluids
   P_crit(k) = CoolProp.Props1SI('Pcrit',fluid{k});
   P_min(k) = 1.05*CoolProp.PropsSI('P','T',T_o,'Q',0,fluid{k});
   P_max(k) = 0.95*P_crit(k);
   
   P(k,:) = P_min(k):(P_max(k)-P_min(k))/(div_P-1):P_max(k);
   
   u_o(k) = CoolProp.PropsSI('U','P',P_o,'T',T_o,fluid{k}); % Estado morto
   v_o(k) = 1/CoolProp.PropsSI('D','P',P_o,'T',T_o,fluid{k});
   s_o(k) = CoolProp.PropsSI('S','P',P_o,'T',T_o,fluid{k});
end

for k = 1:n_fluids
    
    for j=1:div_P % Variação da pressão considerando estado de líquido saturado
        for i = 1:div_Q
            T(k,j,i) = CoolProp.PropsSI('T','P',P(k,j),'Q',x(i),fluid{k})-0.1;   % [K]
            u(k,j,i) = CoolProp.PropsSI('U','P',P(k,j),'Q',x(i),fluid{k});     % [J/kg]
            D(k,j,i) = CoolProp.PropsSI('D','P',P(k,j),'Q',x(i),fluid{k});     % [kg/m3]
            v(k,j,i) = 1./D(k,j,i);                                               % [m3/kg]
            s(k,j,i) = CoolProp.PropsSI('S','P',P(k,j),'Q',x(i),fluid{k});     % [J/kg-K]
            
            ex_u(k,j,i) = ((u(k,j,i) - u_o(k)) ./ 3600) ./(v(k,j)*1000) ;
            ex_v(k,j,i) = ((P_o*(v(k,j,i)-v_o(k))) ./ 3600) ./(v(k,j)*1000) ;
            ex_s(k,j,i) = ((- T_o*(s(k,j,i) - s_o(k))) ./ 3600) ./(v(k,j)*1000) ;
            ex(k,j,i) = u(k,j,i) - u_o(k) + P_o*(v(k,j,i)-v_o(k)) - T_o*(s(k,j,i) - s_o(k)); % Exergia específica [J/kg]
            Ex_kJ_kg(k,j,i) = ex(k,j,i)./1000;              % Conversão para kJ/kg 
            Ex_kJ_l(k,j,i) = (ex(k,j,i)./1000)./(v(k,j)*1000);                % [kJ/L]
            Ex_Wh_kg(k,j,i) = ex(k,j,i)./3600;               % Conversão para Wh/kg 
            Ex_Wh_l(k,j,i) = (ex(k,j,i)./3600)./(v(k,j)*1000);                 % [Wh/L]
            
            for m=1:div_Vol
                C_tank(k,j,i,m) = vessel_cost('Bare Module',P(k,j),Volume(m),'Horizontal','Carbon Steel');
                m_tank(k,j,i,m) = Volume(m)/v(k,j,i);
                C_fluid(k,j,i,m) = fluid_cost(fluid{k},m_tank(k,j,i,m));
                C_st(k,j,i,m) = C_tank(k,j,i,m) + C_fluid(k,j,i,m);
                Ex_st(k,j,i,m) = Ex_Wh_l(k,j,i)*Volume(m);
            end            
        end        
    end
end

u_o_air = CoolProp.PropsSI('U','P',P_o,'T',T_o,'air');
v_o_air = 1/CoolProp.PropsSI('D','P',P_o,'T',T_o,'air');
s_o_air = CoolProp.PropsSI('S','P',P_o,'T',T_o,'air');

%----------------------------------- CAES
P_CAES = P_o:(8000000-P_o)/(div_P-1):8000000;

for j=1:div_P % Variação da pressão considerando estado de líquido saturado
%     T_CAES(j) = T_o;   % [K]
    T_CAES(j) = 330;   % [K]
    u_CAES(j) = CoolProp.PropsSI('U','P',P_CAES(j),'T',T_CAES(j),'air');     % [J/kg]
    D_CAES(j) = CoolProp.PropsSI('D','P',P_CAES(j),'T',T_CAES(j),'air');     % [kg/m3]
    v_CAES(j) = 1/D_CAES(j);                                                    % [m3/kg]
    s_CAES(j) = CoolProp.PropsSI('S','P',P_CAES(j),'T',T_CAES(j),'air');     % [J/kg-K]

    ex_CAES(j) = u_CAES(j) - u_o_air + P_o*(v_CAES(j)-v_o_air) - T_o*(s_CAES(j) - s_o_air); % Exergia específica [J/kg]
    Ex_kJ_kg_CAES(j) = ex_CAES(j)./1000;                                  % Conversão para kJ/kg 
    Ex_kJ_l_CAES(j) = (ex_CAES(j)./1000)/(v_CAES(j)*1000);                % [kJ/L]
    Ex_Wh_kg_CAES(j) = ex_CAES(j)./3600;                                   % Conversão para Wh/kg 
    Ex_Wh_l_CAES(j) = (ex_CAES(j)./3600)/(v_CAES(j)*1000);                 % [Wh/L]
    
    for m=1:div_Vol
        C_st_CAES(j,m) = vessel_cost('Bare Module',P_CAES(j),Volume(m),'Horizontal','Carbon Steel');
        Ex_st_CAES(j,m) = (Ex_Wh_l_CAES(j))*Volume(m); % Exergy stored [kWh] (Wh/l = kWh/m3)
    end
end

%---------------------------------- LAES
% Densidade de exergia para LAES
P_LAES = P_o;
T_LAES = 78;   % [K]
u_LAES = CoolProp.PropsSI('U','P',P_LAES,'T',T_LAES,'air');     % [J/kg]
D_LAES = CoolProp.PropsSI('D','P',P_LAES,'T',T_LAES,'air');     % [kg/m3]
v_LAES = 1/D_LAES;                                              % [m3/kg]
s_LAES = CoolProp.PropsSI('S','P',P_LAES,'T',T_LAES,'air');     % [J/kg-K]

ex_LAES = u_LAES - u_o_air + P_o*(v_LAES-v_o_air) - T_o*(s_LAES - s_o_air); % Exergia específica [J/kg]
Ex_kJ_kg_LAES = ex_LAES./1000;                                  % Conversão para kJ/kg 
Ex_kJ_l_LAES = (ex_LAES./1000)/(v_LAES*1000);                % [kJ/L]
Ex_Wh_kg_LAES = ex_LAES/3600;                                   % Conversão para Wh/kg 
Ex_Wh_l_LAES = (ex_LAES/3600)/(v_LAES*1000);                 % [Wh/L]

%-----------------------------------
P_kPa = P./1000;

%% Figura exergia em função do fluido, pressão e título

for k=1:length(fluid)    % Figura da variação da densidade exergética [Wh/L] com a pressão e título
    figure('color',[1 1 1]);
    hold all
    for i=1:length(x)
        plot(P_kPa(k,:),Ex_Wh_l(k,:,i))
    end
    grid on;
    title(fluid{k})
    xlabel('P [kPa]')
    ylabel('Exergy density [Wh/L]')
    legend(num2str(x'));
    xtickformat('%,1.0f')
    applystyle2plot()
end

for k=1:length(fluid)    % Figura da variação da densidade exergética [Wh/kg] com a pressão e título
    figure('color',[1 1 1]);
    hold all
    for i=1:length(x)
        plot(P_kPa(k,:),Ex_Wh_kg(k,:,i))
    end
    grid on;
    title(fluid{k})
    xlabel('P [kPa]')
    ylabel('Exergy density [Wh/kg]')
    legend(num2str(x'));
    applystyle2plot()
end
% grid on;
% xtickformat('%,1.0f')
% xlabel('P [kPa]')
% ylabel('Exergy density [Wh/L]')
% applystyle2plot()

%% Composição da exergia dos fluidos orgânicos
for k=1:n_fluids    % Figura da variação da densidade exergética [Wh/L] com a pressão e título
    figure('color',[1 1 1]);
    hold all
    plot(P_kPa(k,:),Ex_Wh_l(k,:,1))
    plot(P_kPa(k,:),ex_u(k,:,1))
    plot(P_kPa(k,:),ex_v(k,:,1))
    plot(P_kPa(k,:),ex_s(k,:,1))    
    grid on;
    title(fluid{k})
    xlabel('P [kPa]')
    ylabel('Exergy density [Wh/L]')
    legend('Total specific \chi','\chi_u','\chi_v','\chi_s');
    xtickformat('%,1.0f')
    applystyle2plot()
end

% for k=1:n_fluids % T-s Diagramm
%     plotTS(fluid{k})
%     plot(s(k,:,1)./1000,T(k,:,1))
%         
% end

%% Comparação exergia fluidos orgânicos com CAES e LAES

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
plot(P_CAES./1000,Ex_Wh_l_CAES,'--k')
% plot(P_CAES./1000,ones(length(P_CAES),1)*Ex_Wh_l_LAES,'-.k')
grid on;
xtickformat('%,1.0f')
xlabel('P [kPa]')
ylabel('Exergy density [kWh m^{-3}]')
% legenda1 = cat(2,fluid,{'Compressed Air','Liquid Air'});
legenda1 = cat(2,fluid,{'Compressed Air'});
legend(legenda1,'Location','SouthEast')
% legend(fluid_std,'Location','SouthEast')
applystyle2plot()

figure('color',[1 1 1]);
hold all
for k=1:length(fluid)   % Figura da variação da densidade exergética [Wh/kg] com a pressão
    if k<=length(fluid)/2
        plot(P_kPa(k,:),Ex_Wh_kg(k,:,1)./1000,'Color',k/(length(fluid)/2+1)*[1 1 1])
    else
        plot(P_kPa(k,:),Ex_Wh_kg(k,:,1)./1000,'Color',k/(length(fluid)+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',k/(length(fluid)+1)*[1 1 1],'MarkerIndices',1:10:size(P_kPa,2))
    end
end
plot(P_CAES./1000,Ex_Wh_kg_CAES./1000,'--k')
% plot(P_CAES./1000,ones(length(P_CAES),1)*Ex_Wh_kg_LAES,'-.k')
grid on;
xtickformat('%,1.0f')
xlabel('P [kPa]')
ylabel('Exergy density [kWh kg^{-1}]')
legend(legenda1)
% legend(fluid_std,'Location','SouthEast')
applystyle2plot()

%% Cost analysis

for k = 1:n_fluids % Cost of tank
    figure('color',[1 1 1]);
    hold all
%     for j=1:div_P % Variação da pressão considerando estado de líquido saturado
%     for i = 1:div_Q
    for m=1:div_Vol
        if m<=div_Vol/2
            plot(P_kPa(k,:),C_tank(k,:,1,m),'Color',m/(div_Vol/2+1)*[1 1 1])
        else
            plot(P_kPa(k,:),C_tank(k,:,1,m),'Color',m/(div_Vol+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',m/(div_Vol+1)*[1 1 1],'MarkerIndices',1:10:div_P)
        end    
    end
    title(fluid_std{k})
    legend(num2str(Volume'))
    grid on;
    xtickformat('%,1.0f')
    ytickformat('%,1.1f')
    xlabel('P [kPa]')
    ylabel('Cost of the tank [$]')
    applystyle2plot()
end

for k = 1:n_fluids % Cost of fluid
    figure('color',[1 1 1]);
    hold all
%     for j=1:div_P % Variação da pressão considerando estado de líquido saturado
%     for i = 1:div_Q
    for m=1:div_Vol
        if m<=div_Vol/2
            plot(P_kPa(k,:),C_fluid(k,:,1,m),'Color',m/(div_Vol/2+1)*[1 1 1])
        else
            plot(P_kPa(k,:),C_fluid(k,:,1,m),'Color',m/(div_Vol+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',m/(div_Vol+1)*[1 1 1],'MarkerIndices',1:10:div_P)
        end    
    end
    title(fluid_std{k})
    legend(num2str(Volume'))
    grid on;
    xtickformat('%,1.0f')
    ytickformat('%,1.1f')
    xlabel('P [kPa]')
    ylabel('Cost of fluid [$]')
    applystyle2plot()
end

for k = 1:n_fluids % Cost of storage
    figure('color',[1 1 1]);
    hold all
%     for j=1:div_P % Variação da pressão considerando estado de líquido saturado
%     for i = 1:div_Q
    for m=1:div_Vol
        if m<=div_Vol/2
            plot(P_kPa(k,:),C_st(k,:,1,m),'Color',m/(div_Vol/2+1)*[1 1 1])
        else
            plot(P_kPa(k,:),C_st(k,:,1,m),'Color',m/(div_Vol+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',m/(div_Vol+1)*[1 1 1],'MarkerIndices',1:10:div_P)
        end    
    end
    title(fluid_std{k})
    legend(num2str(Volume'))
    grid on;
    xtickformat('%,1.0f')
    ytickformat('%,1.1f')
    xlabel('P [kPa]')
    ylabel('Cost of storage system [$]')
    ylim([])
    applystyle2plot()
    fig = gcf;
    fig.PaperOrientation = 'landscape';
    fig.PaperSize = [6 4.5];
    fig.PaperPositionMode = 'auto';   
end

for k = 1:n_fluids  % Required mass 
    figure('color',[1 1 1]);
    hold all
%     for j=1:div_P % Variação da pressão considerando estado de líquido saturado
%     for i = 1:div_Q
    for m=1:div_Vol
        if m<=div_Vol/2
            plot(P_kPa(k,:),m_tank(k,:,1,m),'Color',m/(div_Vol/2+1)*[1 1 1])
        else
            plot(P_kPa(k,:),m_tank(k,:,1,m),'Color',m/(div_Vol+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',m/(div_Vol+1)*[1 1 1],'MarkerIndices',1:10:div_P)
        end    
    end
    title(fluid_std{k})
    legend(num2str(Volume'))
    grid on;
    xtickformat('%,1.0f')
    ytickformat('%,1.1f')
    xlabel('P [kPa]')
    ylabel('Mass [kg]')
    applystyle2plot()
end

for k = 1:n_fluids  % Specific volume
    figure('color',[1 1 1]);
    hold all
%     for j=1:div_P % Variação da pressão considerando estado de líquido saturado
%     for i = 1:div_Q
    plot(P_kPa(k,:),v(k,:,1))
    title(fluid_std{k})
    grid on;
    xtickformat('%,1.0f')
    ytickformat('%,1.1f')
    xlabel('P [kPa]')
    ylabel('v [m^3/kg]')
    applystyle2plot()
end

for k = 1:n_fluids % Cost per unit energy
    figure('color',[1 1 1]);
    hold all
%     for j=1:div_P % Variação da pressão considerando estado de líquido saturado
%     for i = 1:div_Q
    for m=1:div_Vol
        if m<=div_Vol/2
            plot(P_kPa(k,:),C_st(k,:,1,m)./Ex_st(k,:,1,m),'Color',m/(div_Vol/2+1)*[1 1 1])
        else
            plot(P_kPa(k,:),C_st(k,:,1,m)./Ex_st(k,:,1,m),'Color',m/(div_Vol+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',m/(div_Vol+1)*[1 1 1],'MarkerIndices',1:10:div_P)
        end    
    end
    title(fluid_std{k})
    legend(num2str(Volume'))
    grid on;
    xtickformat('%,1.0f')
    ytickformat('%,1.1f')
    xlabel('P [kPa]')
    ylabel('Cost per unit energy [$ kWh^{-1}]')
    applystyle2plot()
    fig = gcf;
    fig.PaperOrientation = 'landscape';
    fig.PaperSize = [6 4.5];
    
    ylim([0 10000])
end

%----------------------------------- CAES
figure('color',[1 1 1]);
hold all
for i=1:div_Vol
    if i<=div_Vol/2
        plot(P_CAES./1000,C_st_CAES(:,i),'Color',i/(div_Vol/2+1)*[1 1 1])
    else
        plot(P_CAES./1000,C_st_CAES(:,i),'Color',i/(div_Vol+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',i/(div_Vol+1)*[1 1 1],'MarkerIndices',1:10:div_P)
    end
%     yyaxis left
%     plot(P_CAES./1000,C_st_CAES(:,i))
%     yyaxis right
%     plot(P_CAES./1000,Ex_st(:,i))
end
legend(num2str(Volume'))
grid on;
xtickformat('%,1.0f')
ytickformat('%,1.1f')
xlabel('P [kPa]')
ylabel('Cost of storage [$]')
applystyle2plot()
fig = gcf;
fig.PaperOrientation = 'landscape';
fig.PaperSize = [6 4.5];

figure('color',[1 1 1]);  % Cost per unit energy of CAES
hold all
for i=1:div_Vol
    if i<=div_Vol/2
        plot(P_CAES./1000,C_st_CAES(:,i)./Ex_st_CAES(:,i),'Color',i/(div_Vol/2+1)*[1 1 1])
    else
        plot(P_CAES./1000,C_st_CAES(:,i)./Ex_st_CAES(:,i),'Color',i/(div_Vol+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',i/(div_Vol+1)*[1 1 1],'MarkerIndices',1:10:div_P)
    end
%     yyaxis left
%     plot(P_CAES./1000,C_st_CAES(:,i))
%     yyaxis right
%     plot(P_CAES./1000,Ex_st(:,i))
end
legend(num2str(Volume'))
grid on;
xtickformat('%,1.0f')
ytickformat('%,1.1f')
xlabel('P [kPa]')
ylabel('Cost per unit energy [$ kWh^{-1}]')
ylim([0 10000])
applystyle2plot()
fig = gcf;
fig.PaperOrientation = 'landscape';
fig.PaperSize = [6 4.5];