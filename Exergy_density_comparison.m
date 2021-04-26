clear;
clc;
close all

div_P = 100; % Número de divisões na figura
x=0;
div_Q = length(x);

% FLUIDOS AVALIADOS NO ARTIGO
fluid = {'R152a','R134a','R142b','R365mfc','R236ea','R141b'};
fluid_std = {'R-152a','R-134a','R-142b','R-365mfc','R-236ea','R-141b'};

% Inicialização das variáveis
P_crit = zeros(length(fluid),1);
P_max = zeros(length(fluid),1);

P = zeros(length(fluid),div_P);
T = zeros(length(fluid),div_P);
u = zeros(length(fluid),div_P);
D = zeros(length(fluid),div_P);
v = zeros(length(fluid),div_P);
s = zeros(length(fluid),div_P);

P_CAES = zeros(length(fluid),div_P);
T_CAES = zeros(length(fluid),div_P);
u_CAES = zeros(length(fluid),div_P);
D_CAES = zeros(length(fluid),div_P);
v_CAES = zeros(length(fluid),div_P);
s_CAES = zeros(length(fluid),div_P);

% Estado morto
P_o = 101325;
T_o = 298.15;

% Definição dos limites de pressão e vetor de valores de pressão para cada
% fluido
P_min = 600000;
for k=1:length(fluid)
   P_crit(k) = CoolProp.Props1SI('Pcrit',fluid{k});
   P_max(k) = 0.9*P_crit(k);
   
   P(k,:) = P_min:(P_max(k)-P_min)/(div_P-1):P_max(k);
end

P_CAES = P;

Ex_Wh_l = zeros(length(fluid),size(P,2));
Ex_Wh_kg = zeros(length(fluid),size(P,2));
Ex_Wh_l_CAES = zeros(length(fluid),size(P,2));
Ex_Wh_kg_CAES = zeros(length(fluid),size(P,2));

u_o_air = CoolProp.PropsSI('U','P',P_o,'T',T_o,'air');
v_o_air = 1/CoolProp.PropsSI('D','P',P_o,'T',T_o,'air');
s_o_air = CoolProp.PropsSI('S','P',P_o,'T',T_o,'air');

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

Fig_Comp = figure('color',[1 1 1]);
subplot(2,2,1)
hold on
grid on
grid minor
box on
ylabel({'\chi/\chi_{Comp.Air}';'[(Wh/kg) / (Wh/kg)]'})
xtickformat('%,1.0f')
ytickformat('%.1f')
subplot(2,2,2)
hold on
grid on
grid minor
box on
ylabel({'\chi/\chi_{Comp.Air}';'[(Wh/L) / (Wh/L)]'})
xtickformat('%,1.0f')
subplot(2,2,3)
hold on
grid on
grid minor
box on
ylabel({'\chi/\chi_{Liq.Air}';'[(Wh/kg) / (Wh/kg)]'})
xtickformat('%,1.0f')
ytickformat('%.2f')
subplot(2,2,4)
hold on
grid on
grid minor
box on
ylabel({'\chi/\chi_{Liq.Air}';'[(Wh/L) / (Wh/L)]'})
xtickformat('%,1.0f')
ytickformat('%.2f')

for k = 1:length(fluid)
    u_o(k) = CoolProp.PropsSI('U','P',P_o,'T',T_o,fluid{k}); % Estado morto para o fluido
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
            Ex_kJ_l(k,j,i) = (ex(k,j,i)./1000)/(v(k,j,i)*1000);                % [kJ/L]
            Ex_Wh_kg(k,j,i) = ex(k,j,i)/3600;               % Conversão para Wh/kg 
            Ex_Wh_l(k,j,i) = (ex(k,j,i)/3600)/(v(k,j,i)*1000);                 % [Wh/L]
        end
%         j
        T_CAES(k,j) = T_o;   % [K]
        u_CAES(k,j) = CoolProp.PropsSI('U','P',P_CAES(k,j),'T',T_o,'air');     % [J/kg]
        D_CAES(k,j) = CoolProp.PropsSI('D','P',P_CAES(k,j),'T',T_o,'air');     % [kg/m3]
        v_CAES(k,j) = 1/D_CAES(k,j);                                              % [m3/kg]
        s_CAES(k,j) = CoolProp.PropsSI('S','P',P_CAES(k,j),'T',T_o,'air');     % [J/kg-K]
        
        ex_CAES(k,j) = u_CAES(k,j) - u_o_air + P_o*(v_CAES(k,j)-v_o_air) - T_o*(s_CAES(k,j) - s_o_air);
        Ex_kJ_kg_CAES(k,j) = ex_CAES(k,j)./1000;              % Conversão para kJ/kg 
        Ex_kJ_l_CAES(k,j) = (ex_CAES(k,j)./1000)/(v_CAES(k,j)*1000);                % [kJ/L]
        Ex_Wh_kg_CAES(k,j) = ex_CAES(k,j)./3600;               % Conversão para Wh/kg
        Ex_Wh_l_CAES(k,j) = (ex_CAES(k,j)/3600)/(v_CAES(k,j)*1000);                 % [Wh/L]
    end
    
    subplot(2,2,1)
    if k<=length(fluid)/2
        plot(P(k,:)./1000,Ex_Wh_kg(k,:,1)./Ex_Wh_kg_CAES(k,:),'Color',k/(length(fluid)/2+1)*[1 1 1])
    else
%         plot(P(k,:)./1000,Ex_Wh_kg(k,:,1)./Ex_Wh_kg_CAES(k,:),'--','Color',k/(length(fluid)+1)*[1 1 1])
        plot(P(k,:)./1000,Ex_Wh_kg(k,:,1)./Ex_Wh_kg_CAES(k,:),'Color',k/(length(fluid)+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',k/(length(fluid)+1)*[1 1 1],'MarkerIndices',1:10:size(P,2))
    end
    subplot(2,2,2)
    if k<=length(fluid)/2
        plot(P(k,:)./1000,Ex_Wh_l(k,:,1)./Ex_Wh_l_CAES(k,:),'Color',k/(length(fluid)/2+1)*[1 1 1])
    else
%         plot(P(k,:)./1000,Ex_Wh_l(k,:,1)./Ex_Wh_l_CAES(k,:),'--','Color',k/(length(fluid)+1)*[1 1 1])
        plot(P(k,:)./1000,Ex_Wh_l(k,:,1)./Ex_Wh_l_CAES(k,:),'Color',k/(length(fluid)+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',k/(length(fluid)+1)*[1 1 1],'MarkerIndices',1:10:size(P,2))
    end
    subplot(2,2,3)
    if k<=length(fluid)/2
        plot(P(k,:)./1000,Ex_Wh_kg(k,:,1)./Ex_Wh_kg_LAES,'Color',k/(length(fluid)/2+1)*[1 1 1])
    else
%         plot(P(k,:)./1000,Ex_Wh_kg(k,:,1)./Ex_Wh_kg_LAES,'--','Color',k/(length(fluid)+1)*[1 1 1])
        plot(P(k,:)./1000,Ex_Wh_kg(k,:,1)./Ex_Wh_kg_LAES,'Color',k/(length(fluid)+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',k/(length(fluid)+1)*[1 1 1],'MarkerIndices',1:10:size(P,2))
    end
    
    subplot(2,2,4)
    if k<=length(fluid)/2
        plot(P(k,:)./1000,Ex_Wh_l(k,:,1)./Ex_Wh_l_LAES,'Color',k/(length(fluid)/2+1)*[1 1 1])
    else
%         plot(P(k,:)./1000,Ex_Wh_l(k,:,1)./Ex_Wh_l_LAES,'--','Color',k/(length(fluid)+1)*[1 1 1])
        plot(P(k,:)./1000,Ex_Wh_l(k,:,1)./Ex_Wh_l_LAES,'Color',k/(length(fluid)+1)*[1 1 1],'Marker','o','MarkerSize',2,'MarkerFaceColor',k/(length(fluid)+1)*[1 1 1],'MarkerIndices',1:10:size(P,2))
    end
    
end
xlabel('P [kPa]')
legend(fluid_std)
% Fig_Comp.PaperSize = [6,4.5];