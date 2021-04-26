%Simulação transiente de descarregamento de sistema CAES
%   com base nos dados de Liu e Wang (2016) "A comparative research of two 
%   adiabatic compressed air energy storage systems"
fluid = 'air';

dt=10;

T_amb = 303.15;
P_amb = 101325;

eta_t = 0.9;

P_max = 7200000;
P_min = 4200000;

V_Tank = 1000;

T_Turb = 575;
P_turb_out = 625000;

P_Tank(1) = P_max;
T_Tank(1) = 303.15;
rho_Tank(1) = CoolProp.PropsSI('D','P',P_Tank(1),'T',T_Tank(1),fluid);
u_Tank(1) = CoolProp.PropsSI('U','P',P_Tank(1),'T',T_Tank(1),fluid);
m_Tank(1) = rho_Tank(1)*V_Tank;
U_Tank(1) = m_Tank(1)*u_Tank(1);

%------------------------ DISCHARGING ----------------------------
m_disch_dot = 2.791; % kg/s

T_7(1) = 373;

P_11(1) = P_Tank(1);
T_6(1) = T_Tank(1);
h_11(1) = CoolProp.PropsSI('H','P',P_11(1),'T',T_6(1),fluid);

P_7(1) = P_min;
T_7(1) = T_6(1);
h_7(1) = CoolProp.PropsSI('H','P',P_7(1),'T',T_7(1),fluid);

P_8(1) = P_7(1);
T_8(1) = T_Turb;
h_8(1) = CoolProp.PropsSI('H','P',P_8(1),'T',T_8(1),fluid);
s_8(1) = CoolProp.PropsSI('S','P',P_8(1),'T',T_8(1),fluid);

P_9(1) = P_turb_out;
h_9s(1) = CoolProp.PropsSI('H','P',P_9(1),'S',s_8(1),fluid);
h_9(1) = h_8(1) - eta_t*(h_8(1) - h_9s(1));

P_10(1) = P_9(1);
T_10(1) = T_Turb;
h_10(1) = CoolProp.PropsSI('H','P',P_10(1),'T',T_10(1),fluid);
s_10(1) = CoolProp.PropsSI('S','P',P_10(1),'T',T_10(1),fluid);

P_11(1) = P_amb;
h_11s(1) = CoolProp.PropsSI('H','P',P_11(1),'S',s_10(1),fluid);
h_11(1) = h_10(1) - eta_t*(h_10(1) - h_11s(1));

q_h1(1) = h_8(1) - h_7(1);
q_h2(1) = h_10(1) - h_9(1);

w_t1(1) = h_9(1) - h_8(1);
w_t2(1) = h_11(1) - h_10(1);

t(1) = 0;
i=1;
while P_Tank(i)>P_min
    i=i+1;
    t(i) = t(i-1) + dt;
    m_Tank(i) = m_Tank(i-1) - m_disch_dot*dt;
    rho_Tank(i) = m_Tank(i)/V_Tank;
    U_Tank(i) = U_Tank(i-1) - m_disch_dot*h_11(i-1)*dt;
    u_Tank(i) = U_Tank(i)/m_Tank(i);
    
    T_Tank(i) = T_Tank(i-1);  %CoolProp.PropsSI('T','U',u_T(i),'D',rho_T(i),fluid);
    P_Tank(i) = CoolProp.PropsSI('P','U',u_Tank(i),'D',rho_Tank(i),fluid);
    
    P_11(i) = P_Tank(i);
    T_6(i) = T_Tank(i);
    h_11(i) = CoolProp.PropsSI('H','P',P_11(i),'T',T_6(i),fluid);

    P_7(i) = P_min;
    T_7(i) = T_6(i);
    h_7(i) = CoolProp.PropsSI('H','P',P_7(i),'T',T_7(i),fluid);

    P_8(i) = P_7(i);
    T_8(i) = T_Turb;
    h_8(i) = CoolProp.PropsSI('H','P',P_8(i),'T',T_8(i),fluid);
    s_8(i) = CoolProp.PropsSI('S','P',P_8(i),'T',T_8(i),fluid);

    P_9(i) = P_turb_out;
    h_9s(i) = CoolProp.PropsSI('H','P',P_9(i),'S',s_8(i),fluid);
    h_9(i) = h_8(i) - eta_t*(h_8(i) - h_9s(i));

    P_10(i) = P_9(i);
    T_10(i) = T_Turb;
    h_10(i) = CoolProp.PropsSI('H','P',P_10(i),'T',T_10(i),fluid);
    s_10(i) = CoolProp.PropsSI('S','P',P_10(i),'T',T_10(i),fluid);

    P_11(i) = P_amb;
    h_11s(i) = CoolProp.PropsSI('H','P',P_11(i),'S',s_10(i),fluid);
    h_11(i) = h_10(i) - eta_t*(h_10(i) - h_11s(i));

    q_h1(i) = h_8(i) - h_7(i);
    q_h2(i) = h_10(i) - h_9(i);

    w_t1(i) = h_8(i) - h_9(i);
    w_t2(i) = h_10(i) - h_11(i);
end

E_t = sum(m_disch_dot*(w_t1+w_t2)*dt); % Energy generated in the turbines [J]
% Energy density in [Wh/L]
rho_E = (E_t/3600)/(V_Tank*1000);