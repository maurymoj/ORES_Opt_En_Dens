clear
clc

% FLUIDOS AVALIADOS NA FASE INICIAL
% Fluidos com baixa temperatura crítica
% fluid_l = {'R134a','R245fa','R152a','R236fa','R227ea','R143a'};
% Fluidos com média temperatura crítica
% fluid_m = {'R123','R245ca','butane','n-Pentane'};
% Fluidos com alta temperatura crítica
% fluid_h = {'Benzene','Toluene','MDM','Cyclohexane'};

% Fluidos com melhor desempenho (em termos de densidade de energia)
% fluid = {'R143a','R152a','R134a','butane','n-Pentane','Cyclohexane','Toluene','Benzene'};

fluids = {'R152a','R134a','R142b','R365mfc','R236ea','R141b'}; % ARTIGO
rho_lim = [375,514,0,491,0,0;   % lim superior
           360,446,0,471,0,0];  % lim inferior
k = 1;
fluid = fluids{k};

V=1;
Dt = 1;

P(1) = 2000000;
x(1) = 0.05;

u_o = CoolProp.PropsSI('U','P',P(1),'Q',x(1),fluid);
rho_o = CoolProp.PropsSI('D','P',P(1),'Q',x(1),fluid);
m(1) = rho_o*V;
m_e = 0.001*m(1);

n_steps = m(1)/m_e;

u(1) = u_o;
E(1) = m(1)*u(1);
rho(1) = rho_o;
T(1) = CoolProp.PropsSI('T','U',u(1),'D',rho(1),fluid);
s(1) = CoolProp.PropsSI('S','U',u(1),'D',rho(1),fluid);
h(1) = CoolProp.PropsSI('H','U',u(1),'D',rho(1),fluid);
Q(1) = CoolProp.PropsSI('Q','P',P(1),'D',rho(1),fluid);

h_e(1) = CoolProp.PropsSI('H','P',P(1),'Q',0,fluid);

t(1) = 0;
test = true;
Dt_test = 100*Dt;
% Dt_test = 
% N_steps = floor(0.95*m/m_e)
for i=1:n_steps
    t(i+1) = t(i) + Dt;
    m(i+1) = m(i) - m_e*Dt;
    E(i+1) = E(i) - Dt*m_e*h_e(i);
    rho(i+1) = m(i+1)/V;
    
    if rho(i+1)>rho_lim(1,k)
        u(i+1) = E(i+1)/m(i+1);
        P(i+1) = CoolProp.PropsSI('P','U',u(i+1),'D',rho(i+1),fluid);
        T(i+1) = CoolProp.PropsSI('T','U',u(i+1),'D',rho(i+1),fluid);
        s(i+1) = CoolProp.PropsSI('S','U',u(i+1),'D',rho(i+1),fluid);
        h(i+1) = CoolProp.PropsSI('H','U',u(i+1),'D',rho(i+1),fluid);
        Q(i+1) = CoolProp.PropsSI('Q','U',u(i+1),'D',rho(i+1),fluid);
        h_e(i+1) = CoolProp.PropsSI('H','P',P(i+1),'Q',0,fluid);
    elseif (rho(i+1)<rho_lim(1,k) && test)
        disp(strcat('rho_lim on~',num2str(i)))
        
        t(i+1) = t(i) + Dt_test;
        m(i+1) = m(i) - m_e*Dt_test;
        E(i+1) = E(i) - Dt_test*m_e*h_e(i);
        rho(i+1) = m(i+1)/V;
        
        u(i+1) = E(i+1)/m(i+1);
        P(i+1) = CoolProp.PropsSI('P','U',u(i+1),'D',rho(i+1),fluid);
        T(i+1) = CoolProp.PropsSI('T','U',u(i+1),'D',rho(i+1),fluid);
        s(i+1) = CoolProp.PropsSI('S','U',u(i+1),'D',rho(i+1),fluid);
        h(i+1) = CoolProp.PropsSI('H','U',u(i+1),'D',rho(i+1),fluid);
        Q(i+1) = CoolProp.PropsSI('Q','U',u(i+1),'D',rho(i+1),fluid);
        h_e(i+1) = CoolProp.PropsSI('H','P',P(i+1),'Q',0,fluid);
        test =  0;
    elseif rho(i+1)<rho_lim(1,k) && ~test
        u(i+1) = E(i+1)/m(i+1);
        P(i+1) = CoolProp.PropsSI('P','U',u(i+1),'D',rho(i+1),fluid);
        T(i+1) = CoolProp.PropsSI('T','U',u(i+1),'D',rho(i+1),fluid);
        s(i+1) = CoolProp.PropsSI('S','U',u(i+1),'D',rho(i+1),fluid);
        h(i+1) = CoolProp.PropsSI('H','U',u(i+1),'D',rho(i+1),fluid);
        Q(i+1) = CoolProp.PropsSI('Q','U',u(i+1),'D',rho(i+1),fluid);
        h_e(i+1) = CoolProp.PropsSI('H','P',P(i+1),'Q',0,fluid);
    end
    if Q(i+1)>1
        break;
    end
end

plotTS(fluid)
plot(s./1000,T)

figure('color',[1 1 1])
yyaxis left
plot(t,u)
ylabel('u')
yyaxis right
plot(t,rho)
hold on
lim_s = rho_lim(1,k)*ones(1,length(t));
plot(t,lim_s,'k')
lim_i = rho_lim(2,k)*ones(1,length(t));
plot(t,lim_i,'k')

figure('color',[1 1 1])
yyaxis left
plot(t,P)
ylabel('P')
yyaxis right
plot(t,rho)
hold on
lim_s = rho_lim(1,k)*ones(1,length(t));
plot(t,lim_s,'k')
lim_i = rho_lim(2,k)*ones(1,length(t));
plot(t,lim_i,'k')

figure('color',[1 1 1])
yyaxis left
plot(t,T)
ylabel('T')
yyaxis right
plot(t,rho)
hold on
lim_s = rho_lim(1,k)*ones(1,length(t));
plot(t,lim_s,'k')
lim_i = rho_lim(2,k)*ones(1,length(t));
plot(t,lim_i,'k')

figure('color',[1 1 1])
yyaxis left
plot(t,h)
ylabel('h')
yyaxis right
plot(t,rho)
hold on
lim_s = rho_lim(1,k)*ones(1,length(t));
plot(t,lim_s,'k')
lim_i = rho_lim(2,k)*ones(1,length(t));
plot(t,lim_i,'k')