P_o = 101325;
T_o = 298;

P_min = 1.02*CoolProp.PropsSI('P','T',T_o,'Q',0,'CO2');
P_max = 0.98*CoolProp.Props1SI('CO2','Pcrit');

P = P_min:(P_max - P_min)/9:P_max;

K_H = 1.2:0.2:3;
K_L = K_H;

eta = zeros(length(K_H),length(K_L));
rho = zeros(length(K_H),length(K_L));
t_g = zeros(length(K_H),length(K_L));

% for i=1:length(P)
    for j = 1:length(K_H)
        for k = 1:length(K_L)
            [eta(j,k),rho(j,k),t_g(j,k)] = ORES_tr('fluid','CO2','P_HPT',P_max,'K_V_HPT',K_H(j),'K_V_LPT',K_L(k));
        end
    end
% end

figure('color',[1 1 1])
contour(K_L,K_H,eta,'k','ShowText','on')
hold on;
contour(K_L,K_H,rho,'color',[0.7 0.7 0.7],'ShowText','on')
grid on
contour(K_L,K_H,t_g==3600,[1,1],'r')
xlabel('V_{LPT}')
ylabel('V_{HPT}')
legend('\eta_{RT}','\rho_E')

%%
P_o = 101325;
T_o = 298;

P_min = 1.02*CoolProp.PropsSI('P','T',T_o,'Q',0,'CO2');
P_max = 0.98*CoolProp.Props1SI('CO2','Pcrit');

P = P_min:(P_max - P_min)/9:P_max;

% K_H = 1.2:0.2:3;
% K_L = K_H;
eta = zeros(length(P),1);
rho = zeros(length(P),1);
t_g = zeros(length(P),1);

for i=1:length(P)
    [eta(i),rho(i),t_g(i)] = ORES_tr('fluid','CO2','P_HPT',P(i),'K_V_HPT',12,'K_V_LPT',12);
end

figure('color',[1 1 1])
plot(P,eta)
yyaxis right
plot(P,rho,'r')
legend('\eta_{RT}','\rho_E')
% contour(K_L,K_H,eta,'k','ShowText','on')
% hold on;
% contour(K_L,K_H,rho,'color',[0.7 0.7 0.7],'ShowText','on')
% grid on
% contour(K_L,K_H,t_g==3600,[1,1],'r')
% xlabel('V_{LPT}')
% ylabel('V_{HPT}')
% legend('\eta_{RT}','\rho_E')