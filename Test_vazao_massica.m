clear;
clc;
addpath('C:\Users\Maury-D\Documents\MATLAB\COOLPROP');

% fluid = {'R143a','R152a','R134a','butane','n-Pentane','Cyclohexane','Toluene','Benzene'};
fluid = 'R134a';

% Parametric variables
Quality = 0:0.1:0.4; % OK
Volume = 10.^(-3:2);              % 0.001 m^3 (1 L) < V < 100 m^3 (100.000 L)

P_o = 3e5:5.4e5:3e6;              % 300.000 Pa < P_o < 3 MPa
m_exit = 10.^(-3:-1);          % 0.001 m_o < m_e < 0.1 m_o

% Current problem configuration
Dt = 1;
x(1) = Quality(1);
V = Volume(4);

% TESTE V e P
m_e = m_exit(1);

for k=1:length(m_exit)
    figure('color',[1 1 1])
    hold all;
    grid on;
    for j=length(P_o):-1:3
        P(1) = P_o(j);
        P_min = 0.5*P(1);
        %============================ DISCHARGING ================================%
        u_o = CoolProp.PropsSI('U','P',P(1),'Q',x(1),fluid);
        rho_o = CoolProp.PropsSI('D','P',P(1),'Q',x(1),fluid);
        m_o = rho_o*V;
        
        m_e(1) = m_exit(k)*m_o;
        
        n_steps = floor(0.5*m_o/m_e(1));

        m(1) = m_o;
        u(1) = u_o;
        E(1) = m_o*u(1);
        rho(1) = rho_o;
        T(1) = CoolProp.PropsSI('T','U',u(1),'D',rho(1),fluid);
        s(1) = CoolProp.PropsSI('S','U',u(1),'D',rho(1),fluid);
        h(1) = CoolProp.PropsSI('H','U',u(1),'D',rho(1),fluid);
        Q(1) = CoolProp.PropsSI('Q','P',P(1),'D',rho(1),fluid);

        h_e(1) = CoolProp.PropsSI('H','P',P(1),'Q',0,fluid);

        clear t
        t(1) = 0;
        
        for i=1:n_steps
    %i
            t(i+1) = t(i)+Dt;
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
            elseif rho(i) <= 512
                disp('rho lower than 513 kg/m^3.')
                break;
            end
        end    
        yyaxis left
        plot(t,P./1000)
        j
        yyaxis right
        plot(t,m)
    end
    k
end

%%
P(1) = P_o(end);
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
for i=1:n_steps
    %i
    t(i+1) = t(i)+Dt;
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
    end
    
end