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

fluid = {'R152a','R134a','R142b','R365mfc','R236ea','R141b'}; % TESE

% fluid = 'R152a';P = 675000:25000:3500000;
% fluid = 'R134a'; P = 675000:25000:3500000;
%fluid = 'R152a';P = 675000:25000:3500000;
P=zeros(length(fluid),n_div_P);
DeltaT = 0.01:5:40.01;

eta_I = zeros(length(DeltaT),size(P,2),length(fluid));
eta_II = zeros(length(DeltaT),size(P,2),length(fluid));
w_net = zeros(length(DeltaT),size(P,2),length(fluid));
eta_RT = zeros(length(DeltaT),size(P,2),length(fluid));
    
for k=1:length(fluid)
% k=1
    P_min = CoolProp.PropsSI('P','T',298.15,'Q',0,fluid{k});
    P_max = 0.95*CoolProp.Props1SI(fluid{k},'Pcrit');
    P(k,:) = P_min:(P_max - P_min)/(n_div_P-1):P_max;
    

    for i=1:length(P)
        
        for j=1:length(DeltaT)
            
            T = CoolProp.PropsSI('T','P',P(k,i),'Q',1,fluid{k})+DeltaT(j);
            [eta_RT(j,i,k),eta_I(j,i,k),eta_II(j,i,k),q_in(j,i,k),w_net(j,i,k)]=ORES_st('P_HPT',P(k,i),'T_2',T,'fluid',fluid{k});
        end
    end
end
%% Graphs of eta_RT x P x DT_SH for each of the organic fluids

Sty = {'-','--','-.'}';
n_sty = length(Sty);
Sty = repmat(Sty,1,ceil(length(DeltaT)/n_sty));
Sty2 = Sty';
Sty = Sty2(:);
Col = repmat([0:(0.8-0)/(length(DeltaT)/(n_sty+1)):0.8]',n_sty,3);

for k=1:length(fluid)
    figure('color',[1 1 1])
    hold all
    grid on
    for j=1:length(DeltaT)
       plot(P(k,:)./1000,eta_RT(j,:,k),Sty{j},'Color',Col(j,:))
    end
    axis([500,4400,0,1])
    fig = gcf;
    fig.PaperOrientation = 'landscape';
    fig.PaperSize = [6 4.5];
    xlabel('Pressure [kPa]')
    ylabel('\eta_{RT}')
    legend(num2str(floor(DeltaT')),'Location','northwest')
    applystyle2plot()
    name = strcat('eta_R_x_P_HPT_x_DT_SH_',fluid{k});%,'.pdf');
    xtickformat('%,0.0f')
    ytickformat('%,0.1f')
%     print(name,'-dpdf','-bestfit')%,'-r0')
%     print(name,'-depsc')%,'-r0')
%     saveas(gcf,name)
end

%% Comparative of regions
figure('color',[1 1 1])
hold on
grid on
ylim([0 1])

colors = {'r','g','b','y','m','c'};
cor = 0.1:0.15:1;
for k=1:length(fluid)
    P_comp = [P(k,2:end)./1000,fliplr(P(k,2:end)/1000)];
    eta_comp = [eta_RT(end,2:end,k), fliplr(eta_RT(1,2:end,k))];
%     fill(P_comp,eta_comp,colors{k},'EdgeColor',colors{k},'FaceAlpha',.3)
%     fill(P_comp,eta_comp,cor(k).*[1 1 1],'EdgeColor',cor(k).*[1 1 1],'FaceAlpha',.5)
    fill(P_comp,eta_comp,cor(k).*[1 1 1],'FaceAlpha',.5)
end
ytickformat('%.1f')
xtickformat('%,1.0f')
xlabel('Pressure [kPa]')
ylabel('Round-trip efficiency, \eta_{RT}')
legend(fluid)
%% Comparative graph of eta_RT x P for DT_SH=0 of the organic fluids
figure('color',[1 1 1])
hold all
grid on

for k=1:length(fluid)
    
    plot(P(k,:)./1000,eta_RT(1,:,k)) 
    
    
%     applystyle2plot()
%     name = strcat('eta_R_x_P_HPT_x_DT_SH_',fluid{k});%,'.pdf');
%     
%     print(name,'-dpdf','-bestfit')%,'-r0')
%     saveas(gcf,name)
end

axis([500,4300,0,1])

fig = gcf;
fig.PaperOrientation = 'landscape';
fig.PaperSize = [6 4.5];
xlabel('Pressure [kPa]')
ylabel('\eta_{RT}')
legend(fluid,'Location','northwest')

%%
figure('color',[1 1 1])
hold all
grid on
for i=1:length(DeltaT)
   plot(P./1000,eta_RT(i,:)) 
end
axis([500,3500,0,1])

xlabel('Pressure [kPa]')
ylabel('\eta_{RT}')
legend(num2str(DeltaT'))

%%
figure('color',[1 1 1])
%[CS,H] = contour(P./1000,DeltaT,eta_I,'b');
contour(P./1000,DeltaT,eta_I,'k--');
xlabel('Pressure [kPa]')
ylabel('Degree of superheating [K]')
%clabel(CS,H,'manual',700)
%title('\eta_I')
hold on
%figure('color',[1 1 1])
%[CS,H] = 
contour(P./1000,DeltaT,eta_RT,'Color',[0.7 0.7 0.7]);
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
contour(P./1000,DeltaT,w_net,'k');
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