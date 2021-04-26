fluid = 'R134a';
V_HPT = 0.5%:0.75:15;
V_LPT = 0.47%:1.5/2:15;
%V_HPT = 11.5:0.1:12.5;
%V_LPT = 3:0.1:4;

eta_RT = zeros(length(V_HPT),length(V_LPT));
rho_E_L = zeros(length(V_HPT),length(V_LPT));
eta_I = zeros(length(V_HPT),length(V_LPT));
eta_II = zeros(length(V_HPT),length(V_LPT));

for i=1:length(V_HPT)
    for j=1:length(V_LPT)
        %[V_HPT V_LPT]
        [eta_RT(i,j),rho_E_L(i,j),eta_I(i,j),eta_II(i,j)]=ORES_tr('Delta_t',1800,'V_HPT',V_HPT(i),'V_LPT',V_LPT(j),'Plot',1,'fluid',fluid);
    end
end

%% 
figure('color',[1 1 1]);
[CS,H] = contour(V_HPT,V_LPT,eta_I,'b');
xlabel('Volume at the HPT [m^3]')
ylabel('Volume at the LPT [m^3]')
clabel(CS,H,'manual')
hold on
[CS,H] = contour(V_HPT,V_LPT,eta_II,'r');
xlabel('Volume at the HPT [m^3]')
ylabel('Volume at the LPT [m^3]')
clabel(CS,H,'manual')
legend('\eta_I','\eta_{II}')
figure('color',[1 1 1]);
[CS,H] = contour(V_HPT,V_LPT,eta_RT,'b');
xlabel('Volume at the HPT [m^3]')
ylabel('Volume at the LPT [m^3]')
clabel(CS,H,'manual')
hold on
[CS,H] = contour(V_HPT,V_LPT,rho_E_L,'r');
xlabel('Volume at the HPT [m^3]')
ylabel('Volume at the LPT [m^3]')
clabel(CS,H,'manual')
legend('\eta_{RT}','\rho_{E_L}')

%%
figure('color',[1 1 1]);
hold all
grid on
for j =1:2:length(V_LPT)
    plot(V_HPT,eta_I(:,j))
end
xlabel('V_{HPT}')
ylabel('\eta_I')
legend(num2str(V_LPT(1:2:end)'))

figure('color',[1 1 1]);
hold all
grid on
for i =1:2:length(V_HPT)
    plot(V_LPT,eta_I(i,:))
end
xlabel('V_{LPT}')
ylabel('\eta_I')
legend(num2str(V_HPT(1:2:end)'))