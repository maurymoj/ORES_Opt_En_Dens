%% Calculations
n_div_P = 20;
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

% fluid = 'R152a';P = 675000:25000:3500000;
% fluid = 'R134a'; P = 675000:25000:3500000;
%fluid = 'R152a';P = 675000:25000:3500000;
P_H=zeros(length(fluid),n_div_P);
P_L=zeros(length(fluid),n_div_P);

eta_I = zeros(size(P_L,2),size(P_H,2),length(fluid));
eta_II = zeros(size(P_L,2),size(P_H,2),length(fluid));
w_net = zeros(size(P_L,2),size(P_H,2),length(fluid));
eta_RT = zeros(size(P_L,2),size(P_H,2),length(fluid));
    
for k=1:length(fluid)
% k=1
    P_H_min = 2*CoolProp.PropsSI('P','T',298.15,'Q',0,fluid{k});
    P_H_max = 0.95*CoolProp.Props1SI(fluid{k},'Pcrit');
    P_H(k,:) = P_H_min:(P_H_max - P_H_min)/(n_div_P-1):P_H_max;
    
    P_L_min = CoolProp.PropsSI('P','T',298.15,'Q',0,fluid{k});
    P_L_max = 1.8*CoolProp.PropsSI('P','T',298.15,'Q',0,fluid{k});
    P_L(k,:) = P_L_min:(P_L_max - P_L_min)/(n_div_P-1):P_L_max;
    
    for i=1:length(P_H)
        for j=1:length(P_L)
            T = CoolProp.PropsSI('T','P',P_H(k,i),'Q',1,fluid{k})+0.01;
            [eta_I(j,i,k),eta_II(j,i,k),w_net(j,i,k),eta_RT(j,i,k)]=ORES_st('P_HPT',P_H(k,i),'T_2',T,'P_LPT',P_L(k,j),'fluid',fluid{k});
        end
    end
end
%% Graphs of eta_RT x P x DT_SH for each of the organic fluids

Sty = {'-','--','-.'}';
n_sty = length(Sty);
Sty = repmat(Sty,1,ceil(size(P_L,2)/(n_sty+1)));
Sty2 = Sty';
Sty = Sty2(:);
Col = repmat([0:(0.8-0)/(size(P_L,2)/(n_sty+1)):0.8]',n_sty,3);

for k=1:length(fluid)
%     k
    figure('color',[1 1 1])
    hold all
    grid on
    for j=1:size(P_L,2)
%         j
       plot(P_H(k,:)./1000,eta_RT(j,:,k),Sty{j},'Color',Col(j,:))
    end
    axis([500,4300,0,1])
    set(gcf,'PaperPositionMode','auto')
    set(gcf,'PaperOrientation','landscape')
    xlabel('Pressure [kPa]')
    ylabel('\eta_{RT}')
    legend(num2str(P_L'),'Location','northwest')
    applystyle2plot()
    ylim([0,1])
%     name = strcat('eta_R_x_P_HPT_x_P_LPT',fluid{k});%,'.pdf');
%     
%     print(name,'-dpdf','-bestfit')%,'-r0')
%     print(name,'-depsc')%,'-r0')
%     saveas(gcf,name)
end

%% Comparative graph of eta_RT x P for DT_SH=0 of the organic fluids
figure('color',[1 1 1])
hold all
grid on

for k=1:length(fluid)
    
    plot(P_H(k,:)./1000,eta_RT(1,:,k)) 
    
    
%     applystyle2plot()
%     name = strcat('eta_R_x_P_HPT_x_DT_SH_',fluid{k});%,'.pdf');
%     
%     print(name,'-dpdf','-bestfit')%,'-r0')
%     saveas(gcf,name)
end

axis([500,4300,0,1])

set(gcf,'PaperPositionMode','auto')
set(gcf,'PaperOrientation','landscape')
xlabel('Pressure [kPa]')
ylabel('\eta_{RT}')
legend(fluid,'Location','northwest')

%%
figure('color',[1 1 1])
hold all
grid on
for i=1:length(DeltaT)
   plot(P_H./1000,eta_RT(i,:)) 
end
axis([500,3500,0,1])

xlabel('Pressure [kPa]')
ylabel('\eta_{RT}')
legend(num2str(DeltaT'))

%%
figure('color',[1 1 1])
%[CS,H] = contour(P./1000,DeltaT,eta_I,'b');
contour(P_H./1000,DeltaT,eta_I,'k--');
xlabel('Pressure [kPa]')
ylabel('Degree of superheating [K]')
%clabel(CS,H,'manual',700)
%title('\eta_I')
hold on
%figure('color',[1 1 1])
%[CS,H] = 
contour(P_H./1000,DeltaT,eta_RT,'Color',[0.7 0.7 0.7]);
xlabel('Pressure [kPa]')
ylabel('Degree of superheating [K]')
% clabel(CS,H,'manual',700)
%title('\eta_{RT}')

% figure('color',[1 1 1])
% [CS,H] = 
% contour(P./1000,DeltaT,eta_II,'Color',[0.7 0.7 0.7]);
% xlabel('Pressure [kPa]')
% ylabel('Degree of superheating [K]')
% clabel(CS,H,'manual',700)
% title('\eta_{II}')

%figure('color',[1 1 1])
%[CS,H] = 
contour(P_H./1000,DeltaT,w_net,'k');
xlabel('Pressure [kPa]')
ylabel('Degree of superheating [K]')
% clabel(CS,H,'manual',700);
%title('w_{net}')
legend('\eta_I','\eta_{II}','w_{net}')

%%
% Frentes de pareto
figure('color',[1 1 1])
plot(eta_I(:),eta_II(:),'.')
xlabel('\eta_I')
ylabel('\eta_{II}')
figure('color',[1 1 1])
contour(eta_I,eta_II,w_net,'.')
xlabel('\eta_I')
ylabel('\eta_{II}')
zlabel('w_{net} [J/kg]')

%% Diagrama T-s
fluid = 'R134a';
DeltaT_Ts = 10;
P_Ts = 2300000;
T_Ts = CoolProp.PropsSI('T','P',P_Ts,'Q',1,fluid)+DeltaT_Ts;
[eta_I,eta_II,w_net,eta_RT,q_in]=ORES('P_HPT',P_Ts,'T_2',T_Ts,'fluid',fluid,'Plot',1)
axis([1 2.2 275 400])

fluid = 'R152a';
DeltaT_Ts = 20;
P_Ts = 2300000;
T_Ts = CoolProp.PropsSI('T','P',P_Ts,'Q',1,fluid)+DeltaT_Ts;
[eta_I,eta_II,w_net,eta_RT,q_in]=ORES('P_HPT',P_Ts,'T_2',T_Ts,'fluid',fluid,'Plot',1)
axis([1 2.2 275 400])