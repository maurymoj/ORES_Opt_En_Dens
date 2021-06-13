fluid = 'R141b';

eta_RT = zeros(size(x,1),1);
rho_E_L = zeros(size(x,1),1);
CAPEX = zeros(size(x,1),1);
C_fluid = zeros(size(x,1),1);
C_tank = zeros(size(x,1),1);
C_turb = zeros(size(x,1),1); 
C_pump = zeros(size(x,1),1);
C_aux = zeros(size(x,1),1);

for i=1:size(x,1)
    K_V_H = x(i,1);
    K_V_L = x(i,2);
    P_H = x(i,3);    
    [eta_RT(i),rho_E_L(i),~,~,CAPEX(i),C_fluid(i),C_tank(i),C_turb(i),C_pump(i),C_aux(i)]= ORES_tr('K_V_HPT',K_V_H,'K_V_LPT',K_V_L,'P_HPT',P_H,'fluid',fluid);
end

Obj = [eta_RT,rho_E_L];

groups = ['High density','Balanced','High density','Balanced','Balanced',...
'High density','High efficiency','Balanced','High density','High efficiency',...
'High efficiency','High efficiency','Balanced','Balanced','High density',...
'High efficiency','High efficiency','High efficiency'];
%%

figure('Color',[1 1 1])
plot3(obj(strcmpi(groups,"Balanced"),1),obj(strcmpi(groups,"Balanced"),2),CAPEX(strcmpi(groups,"Balanced"),1),'Or')
hold on
grid on
plot3(obj(strcmpi(groups,"High density"),1),obj(strcmpi(groups,"High density"),2),CAPEX(strcmpi(groups,"High density"),1),'Ob')
plot3(obj(strcmpi(groups,"High efficiency"),1),obj(strcmpi(groups,"High efficiency"),2),CAPEX(strcmpi(groups,"High efficiency"),1),'O','Color',[0.93 0.69 0.13])
legend('Balanced','High density','High efficiency')
xlabel('\eta_{RT}')
ylabel('\rho_E')
zlabel('CAPEX')

figure('Color',[1 1 1])
plot3(Obj(:,1),Obj(:,2),CAPEX,'O')
grid on
xlabel('\eta_{RT}')
ylabel('\rho_E')
zlabel('CAPEX')

figure('Color',[1 1 1])
plot3(x(:,1),x(:,2),x(:,3),'O')
plot3(x(strcmpi(groups,"Balanced"),1),x(strcmpi(groups,"Balanced"),2),x(strcmpi(groups,"Balanced"),3),'Or')
hold on
grid on
plot3(x(strcmpi(groups,"High density"),1),x(strcmpi(groups,"High density"),2),x(strcmpi(groups,"High density"),3),'Ob')
plot3(x(strcmpi(groups,"High efficiency"),1),x(strcmpi(groups,"High efficiency"),2),x(strcmpi(groups,"High efficiency"),3),'O','Color',[0.93 0.69 0.13])
grid on
xlabel('K_{V_H}')
ylabel('K_{V_L}')
zlabel('P_H')