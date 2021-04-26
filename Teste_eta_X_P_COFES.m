fluid = {'R152a','R134a','R142b','R365mfc','R236ea','R141b'}; % ARTIGO

for i=1:length(fluid)
    P_LPT = CoolProp.PropsSI('P','T',298,'Q',0,fluid{i});
    P_min = 1.1*P_LPT;
    P_max = 0.95*CoolProp.Props1SI('Pcrit',fluid{i});
    P(i,:) = P_min:(P_max-P_min)/49:P_max;
    
    for j=1:length(P)
        [eta_RT(i,j),eta_I(i,j),~,q_in(i,j),w_net(i,j)] = COFES_st('fluid',fluid{i},'P_HPT',P(i,j),'P_LPT',P_LPT,'DT_SH',0.01);
    end
    
end

figure('color',[1 1 1])
hold all
grid on
for i=1:length(fluid)
    plot(P(i,:)./1000,eta_I(i,:));
%     [hAx,hLine1,hLine2] = plotyy(P./1000,eta_RT(i,:),P./1000,w_net(i,:));
end
xlabel('P [kPa]')
ylabel('\eta_{I}')
legend(fluid)
% set(groot,'defaultAxesColorOrder',co)
% set(groot,'defaultaxeslinestyleorder',{'-','-','-','--','--','--'})
% axis([0 4000 0 1])

figure('color',[1 1 1])
hold all
grid on
for i=1:length(fluid)
    plot(P(i,:)./1000,eta_RT(i,:));
%     [hAx,hLine1,hLine2] = plotyy(P./1000,eta_RT(i,:),P./1000,w_net(i,:));
end
xlabel('P [kPa]')
ylabel('\eta_{RT}')
legend(fluid)
axis([0 4000 0 1])
% ylabel(hAx(1),'\eta_{RT}') % left y-axis 
% ylabel(hAx(2),'w_{net}') % right y-axis

figure('color',[1 1 1])
hold all
grid on
for i=1:length(fluid)
   plot(P(i,:)./1000,w_net(i,:)./1000);
end
xlabel('P [kPa]')
ylabel('w_{net} [kJ/kg]')
legend(fluid)

figure('color',[1 1 1])
hold all
grid on
for i=1:length(fluid)
   plot(P(i,:)./1000,q_in(i,:)./1000);
end
xlabel('P [kPa]')
ylabel('q_{in} [kJ/kg]')
legend(fluid)