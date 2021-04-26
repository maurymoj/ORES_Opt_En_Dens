% Mass and energy analysis
figure('color',[1 1 1])
subplot(2,2,1)
plot(t_charg,m_HPT(n_steps_disch+1:end),'k'); hold on;
plot(t_charg,m_LPT(n_steps_disch+1:end),'color',[0.7 0.7 0.7]); 
grid on;
ylabel('Total mass [kg]')
legend('m_{HPT}','m_{LPT}','location','east')

subplot(2,2,2)
plot(t_charg,rho_HPT(n_steps_disch+1:end),'k'); hold on; 
plot(t_charg,rho_LPT(n_steps_disch+1:end),'color',[0.7 0.7 0.7]); 
grid on;
ylabel('Density [kg/m^3]')
xlabel('Time [s]')
legend('\rho_{HPT}','\rho_{LPT}','location','east')

subplot(2,2,3)
plot(t_charg,u_HPT(n_steps_disch+1:end),'k'); hold on; 
plot(t_charg,u_LPT(n_steps_disch+1:end),'color',[0.7 0.7 0.7]);
grid on;
ylabel('Specific internal energy [J/kg]')
xlabel('Time [s]')
legend('u_{HPT}','u_{LPT}','location','east')

% Thermodynamic properties
try
    figure('color',[1 1 1])
    subplot(3,2,1)
    plot(t_charg,P_HPT(n_steps_disch+1:end)./1000,'k'); hold on;
    plot(t_charg,P_LPT(n_steps_disch+1:end)./1000,'color',[0.7 0.7 0.7]); 
    grid on;
    legend('P_{HPT}','P_{LPT}','location','east')
    ylabel('Pressure [kPa]')

    subplot(3,2,2)
    plot(t_charg,T_HPT(n_steps_disch+1:end),'k'); hold on;
    plot(t_charg,T_LPT(n_steps_disch+1:end),'color',[0.7 0.7 0.7]); 
    grid on;
    ylabel('Temperature [K]')
    legend('T_{HPT}','T_{LPT}','location','east')

    subplot(3,2,3)
    plot(t_charg,x_HPT(n_steps_disch+1:end),'k'); hold on;
    plot(t_charg,x_LPT(n_steps_disch+1:end),'color',[0.7 0.7 0.7]); 
    grid on;
    ylim([0 1])
    ylabel('Quality [-]')
    legend('x_{HPT}','x_{LPT}','location','northwest')

    subplot(3,2,4)
    plot(t_charg,h_HPT(n_steps_disch+1:end),'k'); hold on;
    plot(t_charg,h_LPT(n_steps_disch+1:end),'color',[0.7 0.7 0.7]); 
    grid on;
    ylabel('Enthalpy [J/kg]')
    legend('h_{HPT}','h_{LPT}','location','northwest')

    subplot(3,2,5)
    plot(t_charg,s_HPT(n_steps_disch+1:end),'k'); hold on;
    plot(t_charg,s_LPT(n_steps_disch+1:end),'color',[0.7 0.7 0.7]); 
    grid on;
    ylabel('Entropy [J/kg K]')
    legend('s_{HPT}','s_{LPT}','location','northwest')
catch
   figure('color',[1 1 1])
    subplot(3,2,1)
    plot(t_charg(1:end-1),P_HPT(n_steps_disch+1:end)./1000,'k'); hold on;
    plot(t_charg(1:end-1),P_LPT(n_steps_disch+1:end)./1000,'color',[0.7 0.7 0.7]); 
    grid on;
    legend('P_{HPT}','P_{LPT}','location','east')
    ylabel('Pressure [kPa]')

    subplot(3,2,2)
    plot(t_charg(1:end-1),T_HPT(n_steps_disch+1:end),'k'); hold on;
    plot(t_charg(1:end-1),T_LPT(n_steps_disch+1:end),'color',[0.7 0.7 0.7]); 
    grid on;
    ylabel('Temperature [K]')
    legend('T_{HPT}','T_{LPT}','location','east')

    subplot(3,2,3)
    plot(t_charg(1:end-1),x_HPT(n_steps_disch+1:end),'k'); hold on;
    plot(t_charg(1:end-1),x_LPT(n_steps_disch+1:end),'color',[0.7 0.7 0.7]); 
    grid on;
    ylim([0 1])
    ylabel('Quality [-]')
    legend('x_{HPT}','x_{LPT}','location','northwest')

    subplot(3,2,4)
    plot(t_charg(1:end-1),h_HPT(n_steps_disch+1:end),'k'); hold on;
    plot(t_charg(1:end-1),h_LPT(n_steps_disch+1:end),'color',[0.7 0.7 0.7]); 
    grid on;
    ylabel('Enthalpy [J/kg]')
    legend('h_{HPT}','h_{LPT}','location','northwest')

    subplot(3,2,5)
    plot(t_charg(1:end-1),s_HPT(n_steps_disch+1:end),'k'); hold on;
    plot(t_charg(1:end-1),s_LPT(n_steps_disch+1:end),'color',[0.7 0.7 0.7]); 
    grid on;
    ylabel('Entropy [J/kg K]')
    legend('s_{HPT}','s_{LPT}','location','northwest') 
end

%%

% Plot the T-s diagramm for the fluid with the defined isobaric lines
%     TSDiag=plotTS(fluid,'Iso_P',[P_HPT(n+1),P_LPT(n+1),P_HPT(end),P_LPT(end)]);
plotTS(fluid,'Iso_P',[P_HPT(n_steps_disch+1),P_LPT(n_steps_disch+1),P_HPT(end),P_LPT(end)]);
plot(s_HPT(n_steps_disch+1:end)./1000,T_HPT(n_steps_disch+1:end),'k','LineWidth',2)
handle_charg(1) = plot(s_LPT(n_steps_disch+1:end)./1000,T_LPT(n_steps_disch+1:end),'k','LineWidth',2);

plotTS(fluid,'Process','Pump_comp',P_5(1),h_5(1),P_6(1),eta_p,'plotColor','k');
handle_charg(2) = plotTS(fluid,'Process','Isobaric',P_6(1),s_6(1),s_7(1),'plotColor','k');

plotTS(fluid,'Process','Pump_comp',P_5(end),h_5(end),P_6(end),eta_p,'plotColor',[0.7 0.7 0.7]);
handle_charg(3) = plotTS(fluid,'Process','Isobaric',P_6(end),s_6(end),s_7(end),'plotColor',[0.7 0.7 0.7]);
grid on;
% T-s diagramm during charging phase
legend(handle_charg,{'State at the tanks','Process diagramm at t=0 s',...
    'Process diagramm at t=t_{charge}'},'location','northwest')
ylim([T_LPT(1)-10 T_2(1)+20]);
%xlim([floor(s_LPT(1)./1000