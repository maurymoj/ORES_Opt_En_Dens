clear
clc
% TESTE DA VARIAÇÃO NO ESTADO NOS TANQUES AO LONGO DO CARREGAMENTO/DESC.

% FLUIDOS AVALIADOS NA FASE INICIAL
% Fluidos com baixa temperatura crítica
% fluid_l = {'R134a','R245fa','R152a','R236fa','R227ea','R143a'};
% Fluidos com média temperatura crítica
% fluid_m = {'R122','R245ca','butane','n-Pentane'};
% Fluidos com alta temperatura crítica
% fluid_h = {'Benzene','Toluene','MDM','Cyclohexane'};

% Fluidos com melhor desempenho (em termos de densidade de energia)
% fluid = {'R143a','R152a','R134a','butane','n-Pentane','Cyclohexane','Toluene','Benzene'};

fluids = {'R152a','R134a','R142b','R365mfc','R236ea','R141b'}; % ARTIGO
rho_lim = [360,446,0,471,0,0;   % lim superior
           375,514,0,491,0,0];  % lim inferior
k = 1;
fluid = fluids{k};

Dt = 1;

% Initial state
V_H=1;
P_H(1) = 3200000;
Q_H(1) = 0.05;

V_L=3;
P_L(1) = 1000000;
Q_L(1) = 0.95;

% Thermodynamic properties

% HPT
u_H_o = CoolProp.PropsSI('U','P',P_H(1),'Q',Q_H(1),fluid);
rho_H_o = CoolProp.PropsSI('D','P',P_H(1),'Q',Q_H(1),fluid);
m_H(1) = rho_H_o*V_H;
m_e = 0.001*m_H(1);

n_steps = m_H(1)/m_e;

u_H(1) = u_H_o;
U_H(1) = m_H(1)*u_H(1);
rho_H(1) = rho_H_o;
T_H(1) = CoolProp.PropsSI('T','U',u_H(1),'D',rho_H(1),fluid);
s_H(1) = CoolProp.PropsSI('S','U',u_H(1),'D',rho_H(1),fluid);
h_H(1) = CoolProp.PropsSI('H','U',u_H(1),'D',rho_H(1),fluid);
Q_H(1) = CoolProp.PropsSI('Q','P',P_H(1),'D',rho_H(1),fluid);

h_e(1) = CoolProp.PropsSI('H','P',P_H(1),'Q',0,fluid);

% LPT
u_L_o = CoolProp.PropsSI('U','P',P_L(1),'Q',Q_L(1),fluid);
rho_L_o = CoolProp.PropsSI('D','P',P_L(1),'Q',Q_L(1),fluid);
m_L(1) = rho_L_o*V_L;
m_i = m_e;

u_L(1) = u_L_o;
U_L(1) = m_L(1)*u_L(1);
rho_L(1) = rho_L_o;
T_L(1) = CoolProp.PropsSI('T','U',u_L(1),'D',rho_L(1),fluid);
s_L(1) = CoolProp.PropsSI('S','U',u_L(1),'D',rho_L(1),fluid);
h_L(1) = CoolProp.PropsSI('H','U',u_L(1),'D',rho_L(1),fluid);
Q_L(1) = CoolProp.PropsSI('Q','P',P_L(1),'D',rho_L(1),fluid);

h_i(1) = CoolProp.PropsSI('H','P',P_L(1),'Q',0,fluid);

t(1) = 0;
test = true;
Dt_test = 100*Dt;
% Dt_test = 
% N_steps = floor(0.95*m/m_e)
for i=1:n_steps
%     i
    t(i+1) = t(i) + Dt;
    m_H(i+1) = m_H(i) - m_e*Dt;
%     if i==283
%         i
%     end
    U_H(i+1) = U_H(i) - Dt*m_e*h_e(i);
    u_H(i+1) = U_H(i+1)/m_H(i+1);
    rho_H(i+1) = m_H(i+1)/V_H;    

    m_L(i+1) = m_L(i) + m_i*Dt;
    U_L(i+1) = U_L(i) + Dt*m_i*h_i(i);
    u_L(i+1) = U_L(i+1)/m_L(i+1);
    rho_L(i+1) = m_L(i+1)/V_L;
    
    
    if (rho_H(i+1)<rho_lim(1,k) || rho_H(i+1)>rho_lim(2,k)) && (rho_L(i+1)<rho_lim(1,k) || rho_L(i+1)>rho_lim(2,k))
%         i
        disp('disch - if 1')
        P_H(i+1) = CoolProp.PropsSI('P','U',u_H(i+1),'D',rho_H(i+1),fluid);
        T_H(i+1) = CoolProp.PropsSI('T','U',u_H(i+1),'D',rho_H(i+1),fluid);
        s_H(i+1) = CoolProp.PropsSI('S','U',u_H(i+1),'D',rho_H(i+1),fluid);
        h_H(i+1) = CoolProp.PropsSI('H','U',u_H(i+1),'D',rho_H(i+1),fluid);
        Q_H(i+1) = CoolProp.PropsSI('Q','U',u_H(i+1),'D',rho_H(i+1),fluid);
        h_e(i+1) = CoolProp.PropsSI('H','P',P_H(i+1),'Q',0,fluid);
        
        P_L(i+1) = CoolProp.PropsSI('P','U',u_L(i+1),'D',rho_L(i+1),fluid);
        T_L(i+1) = CoolProp.PropsSI('T','U',u_L(i+1),'D',rho_L(i+1),fluid);
        s_L(i+1) = CoolProp.PropsSI('S','U',u_L(i+1),'D',rho_L(i+1),fluid);
        h_L(i+1) = CoolProp.PropsSI('H','U',u_L(i+1),'D',rho_L(i+1),fluid);
        Q_L(i+1) = CoolProp.PropsSI('Q','U',u_L(i+1),'D',rho_L(i+1),fluid);
        h_i(i+1) = CoolProp.PropsSI('H','P',P_L(i+1),'Q',0,fluid);
    
    elseif (rho_H(i+1)>rho_lim(1,k) && rho_H(i+1)<rho_lim(2,k)) && (rho_L(i+1)<rho_lim(1,k) || rho_L(i+1)>rho_lim(2,k))
%         i
        disp('disch - if 2')
        if (rho_H(i+1)>rho_lim(1,k) && rho_H(i+1)<rho_lim(2,k)) && (rho_H(i)<rho_lim(1,k) || rho_H(i)>rho_lim(2,k))
            disp('update_HPT')
            u_lim_s_H = u_H(i);
            rho_lim_s_H = rho_H(i);
            P_lim_s_H = CoolProp.PropsSI('P','U',u_lim_s_H,'D',rho_lim_s_H,fluid);
            T_lim_s_H = CoolProp.PropsSI('T','U',u_lim_s_H,'D',rho_lim_s_H,fluid);
            x_lim_s_H = CoolProp.PropsSI('Q','U',u_lim_s_H,'D',rho_lim_s_H,fluid);
            h_lim_s_H = CoolProp.PropsSI('H','U',u_lim_s_H,'D',rho_lim_s_H,fluid);
            s_lim_s_H = CoolProp.PropsSI('S','U',u_lim_s_H,'D',rho_lim_s_H,fluid);

            dm = (rho_H(i) - rho_lim(1,k))*V_H;

            m_lim_i_H = m_H(i) - dm;
            rho_lim_i_H = m_lim_i_H/V_H;
            U_lim_i_H = U_H(i+1) - dm*h_e(i);
            u_lim_i_H = U_lim_i_H/m_lim_i_H;
            P_lim_i_H = CoolProp.PropsSI('P','U',u_lim_i_H,'D',rho_lim_i_H,fluid);
            T_lim_i_H = CoolProp.PropsSI('T','U',u_lim_i_H,'D',rho_lim_i_H,fluid);
            x_lim_i_H = CoolProp.PropsSI('Q','U',u_lim_i_H,'D',rho_lim_i_H,fluid);
            h_lim_i_H = CoolProp.PropsSI('H','U',u_lim_i_H,'D',rho_lim_i_H,fluid);
            s_lim_i_H = CoolProp.PropsSI('S','U',u_lim_i_H,'D',rho_lim_i_H,fluid);
        end
        
        P_H(i+1) = P_lim_i_H + (rho_H(i+1)-rho_lim_i_H)/(rho_lim_s_H-rho_lim_i_H)*(P_lim_s_H - P_lim_i_H);
        T_H(i+1) = T_lim_i_H + (rho_H(i+1)-rho_lim_i_H)/(rho_lim_s_H-rho_lim_i_H)*(T_lim_s_H - T_lim_i_H);
        s_H(i+1) = s_lim_i_H + (rho_H(i+1)-rho_lim_i_H)/(rho_lim_s_H-rho_lim_i_H)*(s_lim_s_H - s_lim_i_H);
        h_H(i+1) = h_lim_i_H + (rho_H(i+1)-rho_lim_i_H)/(rho_lim_s_H-rho_lim_i_H)*(h_lim_s_H - h_lim_i_H);
        Q_H(i+1) = x_lim_i_H + (rho_H(i+1)-rho_lim_i_H)/(rho_lim_s_H-rho_lim_i_H)*(x_lim_s_H - x_lim_i_H);
        h_e(i+1) = CoolProp.PropsSI('H','P',P_H(i+1),'Q',0,fluid);
        
        P_L(i+1) = CoolProp.PropsSI('P','U',u_L(i+1),'D',rho_L(i+1),fluid);
        T_L(i+1) = CoolProp.PropsSI('T','U',u_L(i+1),'D',rho_L(i+1),fluid);
        s_L(i+1) = CoolProp.PropsSI('S','U',u_L(i+1),'D',rho_L(i+1),fluid);
        h_L(i+1) = CoolProp.PropsSI('H','U',u_L(i+1),'D',rho_L(i+1),fluid);
        Q_L(i+1) = CoolProp.PropsSI('Q','U',u_L(i+1),'D',rho_L(i+1),fluid);
        h_i(i+1) = CoolProp.PropsSI('H','P',P_L(i+1),'Q',0,fluid);
        
        
    elseif (rho_H(i+1)<rho_lim(1,k) || rho_H(i+1)>rho_lim(2,k)) && (rho_L(i+1)>rho_lim(1,k) && rho_L(i+1)<rho_lim(2,k))
            disp('disch - if 3')
        
        P_H(i+1) = CoolProp.PropsSI('P','U',u_H(i+1),'D',rho_H(i+1),fluid);
        T_H(i+1) = CoolProp.PropsSI('T','U',u_H(i+1),'D',rho_H(i+1),fluid);
        s_H(i+1) = CoolProp.PropsSI('S','U',u_H(i+1),'D',rho_H(i+1),fluid);
        h_H(i+1) = CoolProp.PropsSI('H','U',u_H(i+1),'D',rho_H(i+1),fluid);
        Q_H(i+1) = CoolProp.PropsSI('Q','U',u_H(i+1),'D',rho_H(i+1),fluid);
        h_e(i+1) = CoolProp.PropsSI('H','P',P_H(i+1),'Q',0,fluid);

        if (rho_L(i+1)>rho_lim(1,k) && rho_L(i+1)<rho_lim(2,k)) && (rho_L(i)<rho_lim(1,k) || rho_L(i)>rho_lim(2,k))
            disp('update_LPT')
            u_lim_i_L = u_L(i);
            rho_lim_i_L = rho_L(i);
            P_lim_i_L = CoolProp.PropsSI('P','U',u_lim_i_L,'D',rho_lim_i_L,fluid);
            T_lim_i_L = CoolProp.PropsSI('T','U',u_lim_i_L,'D',rho_lim_i_L,fluid);
            x_lim_i_L = CoolProp.PropsSI('Q','U',u_lim_i_L,'D',rho_lim_i_L,fluid);
            h_lim_i_L = CoolProp.PropsSI('H','U',u_lim_i_L,'D',rho_lim_i_L,fluid);
            s_lim_i_L = CoolProp.PropsSI('S','U',u_lim_i_L,'D',rho_lim_i_L,fluid);

            dm = (rho_lim(2,k) - rho_L(i))*V_L;

            m_lim_s_L = m_L(i) + dm;
            rho_lim_s_L = m_lim_s_L/V_L;
            U_lim_s_L = U_L(i) + dm*h_i(i);
            u_lim_s_L = U_lim_s_L/m_lim_s_L;
            P_lim_s_L = CoolProp.PropsSI('P','U',u_lim_s_L,'D',rho_lim_s_L,fluid);
            T_lim_s_L = CoolProp.PropsSI('T','U',u_lim_s_L,'D',rho_lim_s_L,fluid);
            x_lim_s_L = CoolProp.PropsSI('Q','U',u_lim_s_L,'D',rho_lim_s_L,fluid);
            h_lim_s_L = CoolProp.PropsSI('H','U',u_lim_s_L,'D',rho_lim_s_L,fluid);
            s_lim_s_L = CoolProp.PropsSI('S','U',u_lim_s_L,'D',rho_lim_s_L,fluid);
        end

        P_L(i+1) = P_lim_i_L + (rho_L(i+1)-rho_lim_i_L)/(rho_lim_s_L-rho_lim_i_L)*(P_lim_s_L - P_lim_i_L);
        T_L(i+1) = T_lim_i_L + (rho_L(i+1)-rho_lim_i_L)/(rho_lim_s_L-rho_lim_i_L)*(T_lim_s_L - T_lim_i_L);
        Q_L(i+1) = x_lim_i_L + (rho_L(i+1)-rho_lim_i_L)/(rho_lim_s_L-rho_lim_i_L)*(x_lim_s_L - x_lim_i_L);
        h_L(i+1) = h_lim_i_L + (rho_L(i+1)-rho_lim_i_L)/(rho_lim_s_L-rho_lim_i_L)*(h_lim_s_L - h_lim_i_L);
        s_L(i+1) = s_lim_i_L + (rho_L(i+1)-rho_lim_i_L)/(rho_lim_s_L-rho_lim_i_L)*(s_lim_s_L - s_lim_i_L);
        h_i(i+1) = CoolProp.PropsSI('H','P',P_L(i+1),'Q',0,fluid);

    elseif (rho_H(i+1)>rho_lim(1,k) && rho_H(i+1)<rho_lim(2,k)) && (rho_L(i+1)>rho_lim(1,k) && rho_L(i+1)<rho_lim(2,k))
             disp('disch - if 4')
        if (rho_H(i+1)>rho_lim(1,k) && rho_H(i+1)<rho_lim(2,k)) && (rho_H(i)<rho_lim(1,k) || rho_H(i)>rho_lim(2,k))
            disp('update_HPT 2')
            u_lim_s_H = u_H(i);
            rho_lim_s_H = rho_H(i);
            P_lim_s_H = CoolProp.PropsSI('P','U',u_lim_s_H,'D',rho_lim_s_H,fluid);
            T_lim_s_H = CoolProp.PropsSI('T','U',u_lim_s_H,'D',rho_lim_s_H,fluid);
            x_lim_s_H = CoolProp.PropsSI('Q','U',u_lim_s_H,'D',rho_lim_s_H,fluid);
            h_lim_s_H = CoolProp.PropsSI('H','U',u_lim_s_H,'D',rho_lim_s_H,fluid);
            s_lim_s_H = CoolProp.PropsSI('S','U',u_lim_s_H,'D',rho_lim_s_H,fluid);

            dm = (rho_H(i) - rho_lim(1,k))*V_H;

            m_lim_i_H = m_H(i) - dm;
            rho_lim_i_H = m_lim_i_H/V_H;
            U_lim_i_H = U_H(i) - dm*h_e(i);
            u_lim_i_H = U_lim_i_H/m_lim_i_H;
            P_lim_i_H = CoolProp.PropsSI('P','U',u_lim_i_H,'D',rho_lim_i_H,fluid);
            T_lim_i_H = CoolProp.PropsSI('T','U',u_lim_i_H,'D',rho_lim_i_H,fluid);
            x_lim_i_H = CoolProp.PropsSI('Q','U',u_lim_i_H,'D',rho_lim_i_H,fluid);
            h_lim_i_H = CoolProp.PropsSI('H','U',u_lim_i_H,'D',rho_lim_i_H,fluid);
            s_lim_i_H = CoolProp.PropsSI('S','U',u_lim_i_H,'D',rho_lim_i_H,fluid);
        end

        P_H(i+1) = P_lim_i_H + (rho_H(i+1)-rho_lim_i_H)/(rho_lim_s_H-rho_lim_i_H)*(P_lim_s_H - P_lim_i_H);
        T_H(i+1) = T_lim_i_H + (rho_H(i+1)-rho_lim_i_H)/(rho_lim_s_H-rho_lim_i_H)*(T_lim_s_H - T_lim_i_H);
        Q_H(i+1) = x_lim_i_H + (rho_H(i+1)-rho_lim_i_H)/(rho_lim_s_H-rho_lim_i_H)*(x_lim_s_H - x_lim_i_H);
        h_H(i+1) = h_lim_i_H + (rho_H(i+1)-rho_lim_i_H)/(rho_lim_s_H-rho_lim_i_H)*(h_lim_s_H - h_lim_i_H);
        s_H(i+1) = s_lim_i_H + (rho_H(i+1)-rho_lim_i_H)/(rho_lim_s_H-rho_lim_i_H)*(s_lim_s_H - s_lim_i_H);
        h_e(i+1) = CoolProp.PropsSI('H','P',P_H(i+1),'Q',0,fluid);

        if (rho_L(i+1)>rho_lim(1,k) && rho_L(i+1)<rho_lim(2,k)) && (rho_L(i)<rho_lim(1,k) || rho_L(i)>rho_lim(2,k))
            disp('update_LPT 2')
            u_lim_i_L = u_L(i);
            rho_lim_i_L = rho_L(i);
            P_lim_i_L = CoolProp.PropsSI('P','U',u_lim_i_L,'D',rho_lim_i_L,fluid);
            T_lim_i_L = CoolProp.PropsSI('T','U',u_lim_i_L,'D',rho_lim_i_L,fluid);
            x_lim_i_L = CoolProp.PropsSI('Q','U',u_lim_i_L,'D',rho_lim_i_L,fluid);
            h_lim_i_L = CoolProp.PropsSI('H','U',u_lim_i_L,'D',rho_lim_i_L,fluid);
            s_lim_i_L = CoolProp.PropsSI('S','U',u_lim_i_L,'D',rho_lim_i_L,fluid);

            dm = (rho_lim(2,k) - rho_L(i))*V_L;

            m_lim_s_L = m_L(i) + dm;
            rho_lim_s_L = m_lim_s_L/V_L;
            U_lim_s_L = U_L(i) + dm*h_i(i);
            u_lim_s_L = U_lim_s_L/m_lim_s_L;
            P_lim_s_L = CoolProp.PropsSI('P','U',u_lim_s_L,'D',rho_lim_s_L,fluid);
            T_lim_s_L = CoolProp.PropsSI('T','U',u_lim_s_L,'D',rho_lim_s_L,fluid);
            x_lim_s_L = CoolProp.PropsSI('Q','U',u_lim_s_L,'D',rho_lim_s_L,fluid);
            h_lim_s_L = CoolProp.PropsSI('H','U',u_lim_s_L,'D',rho_lim_s_L,fluid);
            s_lim_s_L = CoolProp.PropsSI('S','U',u_lim_s_L,'D',rho_lim_s_L,fluid);
        end

        P_L(i+1) = P_lim_i_L + (rho_L(i+1)-rho_lim_i_L)/(rho_lim_s_L-rho_lim_i_L)*(P_lim_s_L - P_lim_i_L);
        T_L(i+1) = T_lim_i_L + (rho_L(i+1)-rho_lim_i_L)/(rho_lim_s_L-rho_lim_i_L)*(T_lim_s_L - T_lim_i_L);
        Q_L(i+1) = x_lim_i_L + (rho_L(i+1)-rho_lim_i_L)/(rho_lim_s_L-rho_lim_i_L)*(x_lim_s_L - x_lim_i_L);
        h_L(i+1) = h_lim_i_L + (rho_L(i+1)-rho_lim_i_L)/(rho_lim_s_L-rho_lim_i_L)*(h_lim_s_L - h_lim_i_L);
        s_L(i+1) = s_lim_i_L + (rho_L(i+1)-rho_lim_i_L)/(rho_lim_s_L-rho_lim_i_L)*(s_lim_s_L - s_lim_i_L);
        h_i(i+1) = CoolProp.PropsSI('H','P',P_L(i+1),'Q',0,fluid);
    else
        disp('WTF')
    end    
       
%     elseif (rho_H(i+1)<rho_lim(1,k) && test)
%         disp(strcat('rho_lim on~',num2str(i)))
%         
%         t(i+1) = t(i) + Dt_test;
%         m_H(i+1) = m_H(i) - m_e*Dt_test;
%         E_H(i+1) = E_H(i) - Dt_test*m_e*h_e(i);
%         rho_H(i+1) = m_H(i+1)/V_H;
%         
%         u_H(i+1) = E_H(i+1)/m_H(i+1);
%         P_H(i+1) = CoolProp.PropsSI('P','U',u_H(i+1),'D',rho_H(i+1),fluid);
%         T_H(i+1) = CoolProp.PropsSI('T','U',u_H(i+1),'D',rho_H(i+1),fluid);
%         s_H(i+1) = CoolProp.PropsSI('S','U',u_H(i+1),'D',rho_H(i+1),fluid);
%         h_H(i+1) = CoolProp.PropsSI('H','U',u_H(i+1),'D',rho_H(i+1),fluid);
%         Q_H(i+1) = CoolProp.PropsSI('Q','U',u_H(i+1),'D',rho_H(i+1),fluid);
%         h_e(i+1) = CoolProp.PropsSI('H','P',P_H(i+1),'Q',0,fluid);
%         test =  0;
%     elseif rho_H(i+1)<rho_lim(1,k) && ~test
%         u_H(i+1) = E_H(i+1)/m_H(i+1);
%         P_H(i+1) = CoolProp.PropsSI('P','U',u_H(i+1),'D',rho_H(i+1),fluid);
%         T_H(i+1) = CoolProp.PropsSI('T','U',u_H(i+1),'D',rho_H(i+1),fluid);
%         s_H(i+1) = CoolProp.PropsSI('S','U',u_H(i+1),'D',rho_H(i+1),fluid);
%         h_H(i+1) = CoolProp.PropsSI('H','U',u_H(i+1),'D',rho_H(i+1),fluid);
%         Q_H(i+1) = CoolProp.PropsSI('Q','U',u_H(i+1),'D',rho_H(i+1),fluid);
%         h_e(i+1) = CoolProp.PropsSI('H','P',P_H(i+1),'Q',0,fluid);
%     end
    if Q_H(i)>1
        disp('Q_H > 1')
        break;
    elseif Q_L(i)>1
        disp('Q_L > 1')
        break;
    end
end

[m_H(1) m_L(1);m_H(end) m_L(end)]

[rho_H(1) rho_L(1);rho_H(end) rho_L(end)]

plotTS(fluid)
plot(s_H./1000,T_H)

%%
plotTS(fluid)
plot(s_H./1000,T_H)

figure('color',[1 1 1])
yyaxis left
plot(t,u_H)
ylabel('u')
yyaxis right
plot(t,rho_H)
hold on
lim_s = rho_lim(1,k)*ones(1,length(t));
plot(t,lim_s,'k')
lim_i = rho_lim(2,k)*ones(1,length(t));
plot(t,lim_i,'k')

figure('color',[1 1 1])
yyaxis left
plot(t,P_H)
ylabel('P')
yyaxis right
plot(t,rho_H)
hold on
lim_s = rho_lim(1,k)*ones(1,length(t));
plot(t,lim_s,'k')
lim_i = rho_lim(2,k)*ones(1,length(t));
plot(t,lim_i,'k')

figure('color',[1 1 1])
yyaxis left
plot(t,T_H)
ylabel('T')
yyaxis right
plot(t,rho_H)
hold on
lim_s = rho_lim(1,k)*ones(1,length(t));
plot(t,lim_s,'k')
lim_i = rho_lim(2,k)*ones(1,length(t));
plot(t,lim_i,'k')

figure('color',[1 1 1])
yyaxis left
plot(t,h_H)
ylabel('h')
yyaxis right
plot(t,rho_H)
hold on
lim_s = rho_lim(1,k)*ones(1,length(t));
plot(t,lim_s,'k')
lim_i = rho_lim(2,k)*ones(1,length(t));
plot(t,lim_i,'k')

%%
plotTS(fluid)
plot(s_L./1000,T_L)

figure('color',[1 1 1])
yyaxis left
plot(t,u_L)
ylabel('u')
yyaxis right
plot(t,rho_L)
hold on
lim_s = rho_lim(1,k)*ones(1,length(t));
plot(t,lim_s,'k')
lim_i = rho_lim(2,k)*ones(1,length(t));
plot(t,lim_i,'k')

figure('color',[1 1 1])
yyaxis left
plot(t,P_L)
ylabel('P')
yyaxis right
plot(t,rho_L)
hold on
lim_s = rho_lim(1,k)*ones(1,length(t));
plot(t,lim_s,'k')
lim_i = rho_lim(2,k)*ones(1,length(t));
plot(t,lim_i,'k')

figure('color',[1 1 1])
yyaxis left
plot(t,T_L)
ylabel('T')
yyaxis right
plot(t,rho_L)
hold on
lim_s = rho_lim(1,k)*ones(1,length(t));
plot(t,lim_s,'k')
lim_i = rho_lim(2,k)*ones(1,length(t));
plot(t,lim_i,'k')

figure('color',[1 1 1])
yyaxis left
plot(t,h_L)
ylabel('h')
yyaxis right
plot(t,rho_L)
hold on
lim_s = rho_lim(1,k)*ones(1,length(t));
plot(t,lim_s,'k')
lim_i = rho_lim(2,k)*ones(1,length(t));
plot(t,lim_i,'k')