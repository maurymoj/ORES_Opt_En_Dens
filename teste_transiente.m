figure('color',[1 1 1])
subplot(3,2,1)
%yyaxis left
plot(m_disch)
ylabel('Mass flow rate [kg/s]')
grid on;
%yyaxis right
subplot(3,2,2)
plot(m_HPT(1:end-i))
hold on
plot(m_LPT(1:end-i))
ylabel('m [kg]')
legend('m_{HPT}','m_{LPT}','location','west')
grid on;

subplot(3,2,3)
plot(P_HPT(1:end-i)./1000)
hold on
plot(P_LPT(1:end-i)./1000)
ylabel('P [kPa]')
legend('P_{HPT}','P_{LPT}','location','west')
grid on;

subplot(3,2,4)
plot(T_HPT(1:end-i))
hold on
plot(T_LPT(1:end-i))
ylabel('T [K]')
legend('T_{HPT}','T_{LPT}','location','west')
grid on;

subplot(3,2,5)
%yyaxis left
plot(x_HPT(1:end-i))
hold on
plot(x_LPT(1:end-i))
ylabel('Quality [-]')
ylim([0 1])
legend('x_{HPT}','x_{LPT}','location','northeast')
grid on;
% yyaxis right
subplot(3,2,6)
plot(rho_HPT(1:end-i))
hold on
plot(rho_LPT(1:end-i))
% plot(ones(1,length(h_1))*rho_lim_i_HPT,'g')
% plot(ones(1,length(h_1))*rho_lim_s_HPT,'g')
% plot(ones(1,length(h_1))*rho_lim_i_LPT,'k-')
% plot(ones(1,length(h_1))*rho_lim_s_LPT,'k-')
ylabel('Density [kg/m^3]')
xlabel('Time [s]')
legend('\rho_{HPT}','\rho_{LPT}','location','west')
grid on;

%% ------------------------------------------------------------------------
figure('color',[1 1 1])
subplot(3,2,1)
% yyaxis left
plot(m_charg)
ylabel('mass flow rate [kg/s]')
grid on;
% yyaxis right
subplot(3,2,2)
plot(m_HPT(end-i:end))
hold on
plot(m_LPT(end-i:end))
ylabel('m [kg]')
legend('m_{HPT}','m_{LPT}','location','east')
grid on;

subplot(3,2,3)
% yyaxis left
plot(P_HPT(end-i:end)./1000)
hold on
plot(P_LPT(end-i:end)./1000)
ylabel('P [kPa]')
legend('P_{HPT}','P_{LPT}','location','east')
grid on;
% yyaxis right
subplot(3,2,4)
plot(T_HPT(end-i:end))
hold on
plot(T_LPT(end-i:end))
ylabel('T [K]')
legend('T_{HPT}','T_{LPT}','location','east')
grid on;

subplot(3,2,5)
% yyaxis left
plot(x_HPT(end-i:end))
hold on
plot(x_LPT(end-i:end))
ylabel('Quality [-]')
ylim([0 1])
legend('x_{HPT}','x_{LPT}','location','northwest')
grid on

% yyaxis right
subplot(3,2,6)
plot(rho_HPT(end-i:end))
hold on
plot(rho_LPT(end-i:end))
% plot(ones(1,length(h_7))*rho_lim_i_HPT,'g')
% plot(ones(1,length(h_7))*rho_lim_s_HPT,'g')
% plot(ones(1,length(h_7))*rho_lim_i_LPT,'k-')
% plot(ones(1,length(h_7))*rho_lim_s_LPT,'k-')
ylabel('Density [kg/m^3]')
legend('\rho_{HPT}','\rho_{LPT}','location','east')
grid on;