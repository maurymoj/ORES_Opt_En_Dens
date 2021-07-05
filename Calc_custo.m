fluid = 'R141b';
% fluid = 'R365mfc';
eta_RT = zeros(size(x,1),1);
rho_E_L = zeros(size(x,1),1);
CAPEX = zeros(size(x,1),1);
C_fluid = zeros(size(x,1),1);
C_tank = zeros(size(x,1),1);
C_turb = zeros(size(x,1),1);
C_pump = zeros(size(x,1),1);
C_aux = zeros(size(x,1),1);

for i = 1:size(x,1)
    K_V_H = x(i,1);
    K_V_L = x(i,2);
    P_H = x(i,3);
        
    [eta_RT(i),rho_E_L(i),~,~,CAPEX(i),C_fluid(i),C_tank(i),C_turb(i),C_pump(i),C_aux(i)]= ORES_tr('K_V_HPT',K_V_H,'K_V_LPT',K_V_L,'P_HPT',P_H,'fluid',fluid);
end

Obj = [eta_RT,rho_E_L];

Costs = [C_fluid,C_tank,C_turb,C_pump,C_aux]./1000000;

%%
figure('Color',[1 1 1])
plot3(Obj(:,1),Obj(:,2),CAPEX,'O')
grid on
xlabel('\eta_{RT}')
ylabel('\rho_E')
zlabel('CAPEX')

figure('Color',[1 1 1])
plot3(x(:,1),x(:,2),x(:,3),'O')
grid on
xlabel('K_{V_H}')
ylabel('K_{V_L}')
zlabel('P_H')


lim_RT = prctile(eta_RT,[33 66]);
for i=1:size(x,1)
    if eta_RT(i) < lim_RT(1)
        groups(i) = "High Density";
    elseif eta_RT(i) > lim_RT(2)
        groups(i) = "High efficiency";
    else
        groups(i) = "Balanced";
    end    
end

figure('Color',[1 1 1]) % 2d Figure with objective functions
plot(Obj(strcmpi(groups,"Balanced"),1),Obj(strcmpi(groups,"Balanced"),2),'Or')
hold on
grid on
plot(Obj(strcmpi(groups,"High density"),1),Obj(strcmpi(groups,"High density"),2),'Ob')
plot(Obj(strcmpi(groups,"High efficiency"),1),Obj(strcmpi(groups,"High efficiency"),2),'O','Color',[0.93 0.69 0.13])
legend('Balanced','High density','High efficiency')
xlabel('Round-trip efficiency, \eta_{RT}')
ylabel('Energy density, \rho_{E_v}')

figure('Color',[1 1 1]) % 3d Figure with objective functions and CAPEX
plot3(Obj(strcmpi(groups,"Balanced"),1),Obj(strcmpi(groups,"Balanced"),2),CAPEX(strcmpi(groups,"Balanced"),1),'Or')
hold on
grid on
plot3(Obj(strcmpi(groups,"High density"),1),Obj(strcmpi(groups,"High density"),2),CAPEX(strcmpi(groups,"High density"),1),'Ob')
plot3(Obj(strcmpi(groups,"High efficiency"),1),Obj(strcmpi(groups,"High efficiency"),2),CAPEX(strcmpi(groups,"High efficiency"),1),'O','Color',[0.93 0.69 0.13])
legend('Balanced','High density','High efficiency')
xlabel('\eta_{RT}')
ylabel('\rho_E')
zlabel('CAPEX')

figure('Color',[1 1 1])
plot3(x(strcmpi(groups,"Balanced"),1),x(strcmpi(groups,"Balanced"),2),x(strcmpi(groups,"Balanced"),3),'Or')
hold on
grid on
plot3(x(strcmpi(groups,"High density"),1),x(strcmpi(groups,"High density"),2),x(strcmpi(groups,"High density"),3),'Ob')
plot3(x(strcmpi(groups,"High efficiency"),1),x(strcmpi(groups,"High efficiency"),2),x(strcmpi(groups,"High efficiency"),3),'O','Color',[0.93 0.69 0.13])
xlabel('K_{V_H}')
ylabel('K_{V_L}')
zlabel('P_H')

%% Plot Obj com text
figure('Color',[1 1 1])
plot(Obj(:,1),Obj(:,2),'O')
hold on
grid on
xlim([0.6 0.74])
ylim([0.4 2])
text(Obj(I,1),Obj(I,2) + 0.05,num2str([1:18]'))
xtickformat('%.2f')
ytickformat('%.2f')


xlabel('Round-trip efficiency, \eta_{RT}')
ylabel('Energy density, \rho_{E_v}')

%% Plot CAPEX + Obj + Vars de decisao
figure('Color',[1 1 1])
subplot(3,1,1)
% stacked bar graph for costs
[~,I] = sort(eta_RT(:));
Costs_ord = Costs(I,:);

newcolors = [0.15 0.23 0.37
             0.20 0.30 0.49
             0.31 0.40 0.58
             0.73 0.83 0.96
             0.87 0.92 0.98];
         
% colororder(newcolors)
b = bar(Costs_ord,'stacked','DisplayName','Custos');
legend('C_{fluid}','C_{tank}','C_{Ex}','C_{pump}','C_{aux}','Location','northwestoutside')
grid on
ylabel('Cost (M$)')
ytickformat('%.2f')
ylim([0 8])
xlim([1 18])

for i=1:size(newcolors,1)
    b(i).FaceColor = newcolors(i,:);
    b(i).EdgeColor = newcolors(i,:);
end

subplot(3,1,2)
i=[1:size(Obj,1)]';
[hAx,hLine1,hLine2] = plotyy(i,Obj(I,1),i,Obj(I,2));
hold on
ylabel(hAx(1),'\eta_{RT}') % left y-axis 
ylabel(hAx(2),'\rho_{E_V}, (kWh m^{-3})') % right y-axis
ylim(hAx(1),[0.60 0.76])
ylim(hAx(2),[0 2])
xlim(hAx(1),[1 18])
xlim(hAx(2),[1 18])
ytickformat(hAx(1),'%.2f')
ytickformat(hAx(2),'%.1f')
grid on
legend('\eta_{RT}','\rho_{E_V}','Location','northwestoutside')

subplot(3,1,3)
[hAx,hLine1,hLine2] = plotyy([i i],[x_MPa(I,1) x_MPa(I,2)],i,x_MPa(I,3));
hold on
ylabel(hAx(1),'K_V') % left y-axis 
ytickformat(hAx(1),'%.2f')
ytickformat(hAx(2),'%.1f')
ylabel(hAx(2),'P_H, (MPa)') % right y-axis
ylim(hAx(2),[2.8 4])
xlim(hAx(1),[1 18])
xlim(hAx(2),[1 18])
grid on
legend('K_{V_H}','K_{V_L}','P_H','Location','northwestoutside')