%% Teste variação da densidade em função do título e pressão
fluid = {'R143a','R152a','R134a','butane','n-Pentane','Cyclohexane','Toluene','Benzene'};
% fluid = 'R134a';

P=100000:100000:3000000;
x = 0:0.01:0.1;

for k=1:length(fluid)
    rho = zeros(size(P,1),size(x,1));

    for i=1:length(P)
        for j=1:size(x,2)
            rho(i,j) = CoolProp.PropsSI('D','P',P(i),'Q',x(j),fluid{k});
        end
    end

    figure('color',[1 1 1])
    hold all;
    grid on;

    for j=1:size(x,2)
        plot(P(:)./1000,rho(:,j))
    %     pause(0.1)
    end
    title(fluid{k})
    legend(num2str(x'))
    xlabel('P [kPa]')
    ylabel('\rho [kg/m^3]')
end

%% Teste Descarregamento tanque

clear;
clc;
addpath('C:\Users\Maury-D\Documents\MATLAB\COOLPROP');

% fluid = {'R143a','R152a','R134a','butane','n-Pentane','Cyclohexane','Toluene','Benzene'};
fluids = {'R152a','R134a','R142b','R365mfc','R236ea','R141b'}; % ARTIGO

% Parametric variables
Volume = 10.^(-3:2);              % 1 kg < m_o < 10.000 kg
P_o = 100000:100000:3000000;      % 100 kPa < P_o < 3 MPa
m_exit = 10.^(-3:-1);
Quality = 0:0.1:0.4;
% fluid = 'R134a';
% fluid = 'R152a';
fluid = fluids{1};
% CoolProp.PropsSI('D','P',P_o(1),'Q',0,'R134a')

% Current problem configuration
Dt = 1;
V = Volume(1);
P(1) = P_o(end);
x(1) = Quality(1);
P_min = 0.5*P(1);

%============================ DISCHARGING ================================%
u_o = CoolProp.PropsSI('U','P',P(1),'Q',x(1),fluid);
rho_o = CoolProp.PropsSI('D','P',P(1),'Q',x(1),fluid);
m_o = rho_o*V;
m_e(1) = m_exit(1)*m_o;
n_steps = m_o/m_e(1);

m(1) = m_o;
u(1) = u_o;
E(1) = m_o*u(1);
rho(1) = rho_o;
T(1) = CoolProp.PropsSI('T','U',u(1),'D',rho(1),fluid);
s(1) = CoolProp.PropsSI('S','U',u(1),'D',rho(1),fluid);
h(1) = CoolProp.PropsSI('H','U',u(1),'D',rho(1),fluid);
Q(1) = CoolProp.PropsSI('Q','P',P(1),'D',rho(1),fluid);

h_e(1) = CoolProp.PropsSI('H','P',P(1),'Q',0,fluid);

t(1) = 0;
test = 0;

for i=1:n_steps
%     i
    if i==442
        i
    end
    t(i+1) = t(i) + Dt;
    m_e(i+1) = m_e(i);                   % EXIT MASS FLOW RATE CONSTANT
    m(i+1) = m(i) - Dt*m_e(i);
    E(i+1) = E(i) - Dt*m_e(i)*h_e(i);
    %u(i+1) = u(i) - Dt*( m_e(i)/m(i) )*( u(i) - h_e(i) );
    
    rho(i+1) = m(i+1)/V;
    
    u(i+1) = E(i+1)/m(i+1);
    P(i+1) = CoolProp.PropsSI('P','U',u(i+1),'D',rho(i+1),fluid);
    T(i+1) = CoolProp.PropsSI('T','U',u(i+1),'D',rho(i+1),fluid);
    s(i+1) = CoolProp.PropsSI('S','U',u(i+1),'D',rho(i+1),fluid);
    h(i+1) = CoolProp.PropsSI('H','U',u(i+1),'D',rho(i+1),fluid);
    Q(i+1) = CoolProp.PropsSI('Q','U',u(i+1),'D',rho(i+1),fluid);
    h_e(i+1) = CoolProp.PropsSI('H','P',P(i+1),'Q',0,fluid);    
    
    if P(i+1)<=P_min
        disp('P lower or equal to P_min.')
        break;
    elseif T(i+1) < 293
        disp('T lower than ambient temperature.')
        break;
    elseif test==0 && rho(i+1)<514 && strcmp(fluid,'R134a')
        i
        m(i+1) = m(i) - 0.2*m(i);
        E(i+1) = E(i) - 0.2*m(i)*h_e(i);
%         t(i+1) = t(i) + 10*Dt;
%         m(i+1) = m(i) - 10*Dt*m_e(i);
%         E(i+1) = E(i) - 10*Dt*m_e(i)*h_e(i);
        u(i+1) = E(i+1)/m(i+1);
        rho(i+1) = m(i+1)/V;
        P(i+1) = CoolProp.PropsSI('P','U',u(i+1),'D',rho(i+1),fluid);
        T(i+1) = CoolProp.PropsSI('T','U',u(i+1),'D',rho(i+1),fluid);
        s(i+1) = CoolProp.PropsSI('S','U',u(i+1),'D',rho(i+1),fluid);
        h(i+1) = CoolProp.PropsSI('H','U',u(i+1),'D',rho(i+1),fluid);
        Q(i+1) = CoolProp.PropsSI('Q','U',u(i+1),'D',rho(i+1),fluid);
        h_e(i+1) = CoolProp.PropsSI('H','P',P(i+1),'Q',0,fluid);    
        
        test=1;
    elseif Q(i+1) > 1
        disp('x>1')
        break;
    end
    
end

% taxas médias de variação
dm_dt = (m(end)-m(1))/t(end)
dE_dt = (E(end)-E(1))/t(end)
drho_dt = (rho(end)-rho(1))/t(end)
dhe_dt = (h_e(end)-h_e(1))/t(end)
dP_dt = (P(end)-P(1))/t(end)
%%
figure('color',[1 1 1])
hAx = plotyy(t,m,t,E./1000);
grid on;
ylabel(hAx(1),'m [kg]') % left y-axis 
ylabel(hAx(2),'E [kJ]') % right y-axis
% legend('m','E')

figure('color',[1 1 1])
hAx = plotyy(t,rho,t,u./1000);
grid on;
ylabel(hAx(1),'\rho [kg/m^3]') % left y-axis 
ylabel(hAx(2),'u [kJ/kg]') % right y-axis
% legend('\rho','u')

%% 
figure('color',[1 1 1])
hAx = plotyy(t,Q,t,h_e./1000);
grid on;
ylabel(hAx(1),'x') % left y-axis 
ylabel(hAx(2),'h_e [kJ/kg]') % right y-axis
% legend('Q','h_e')

figure('color',[1 1 1])
hAx = plotyy(t,T,t,P./1000);
grid on;
ylabel(hAx(1),'T [K]') % left y-axis 
ylabel(hAx(2),'P [kPa]') % right y-axis
% legend('T','P')

%% 
figure('color',[1 1 1])
hAx = plotyy(t(1:end-1),Q,t(1:end-1),h_e./1000);
grid on;
ylabel(hAx(1),'x') % left y-axis 
ylabel(hAx(2),'h_e [kJ/kg]') % right y-axis
% legend('Q','h_e')

figure('color',[1 1 1])
hAx = plotyy(t(1:end-1),T,t(1:end-1),P./1000);
grid on;
ylabel(hAx(1),'T [K]') % left y-axis 
ylabel(hAx(2),'P [kPa]') % right y-axis
% legend('T','P')