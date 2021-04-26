fluid = 'R141b';
n_subd = 100;

%% Teste bomba
P1 = 800000;
h1 = CoolProp.PropsSI('H','P',800000,'T',293.15,'R141b');
P2 = 3200000;
eta_p = 0.9;

s_pump = zeros(1,n_subd);
T_pump = zeros(1,n_subd);
h_pump = zeros(1,n_subd);
P_pump = P1:(P2-P1)/(n_subd-1):P2;
h_pump(1) = h1;
% if phase is saturated liquid A elseif subcooled B
s_pump(1) = CoolProp.PropsSI('S','P',P_pump(1),'H',h_pump(1),fluid);
T_pump(1) = CoolProp.PropsSI('T','P',P_pump(1),'H',h_pump(1),fluid);

for j=2:n_subd
    h_s = CoolProp.PropsSI('H','P',P_pump(j),'S',s_pump(j-1),fluid);
    h_pump(j) = pump('h_in',h_pump(j-1),'h_s',h_s,'eta_p',eta_p);
    T_pump(j) = CoolProp.PropsSI('T','P',P_pump(j),'H',h_pump(j),fluid);
    s_pump(j) = CoolProp.PropsSI('S','P',P_pump(j),'H',h_pump(j),fluid);
end

h_s = CoolProp.PropsSI('H','P',P2,'S',s_pump(1),fluid);
h2 = pump('h_in',h1,'h_s',h_s,'eta_p',eta_p);
T2 = CoolProp.PropsSI('T','P',P2,'H',h2,fluid);
s2 = CoolProp.PropsSI('S','P',P2,'H',h2,fluid);


figure('Color',[1 1 1])
plot(s_pump./1000,T_pump)
hold on
plot(s_pump(1)/1000,T_pump(1),'ro')
plot(s2/1000,T2,'ro')

%% Teste bomba
P1 = 3200000;
h1 = CoolProp.PropsSI('H','P',3200000,'T',293.15,'R141b');
P2 = 800000;
eta_p = 0.9;

s_pump = zeros(1,n_subd);
T_pump = zeros(1,n_subd);
h_pump = zeros(1,n_subd);
P_pump = P1:(P2-P1)/(n_subd-1):P2;
h_pump(1) = h1;
% if phase is saturated liquid A elseif subcooled B
s_pump(1) = CoolProp.PropsSI('S','P',P_pump(1),'H',h_pump(1),fluid);
T_pump(1) = CoolProp.PropsSI('T','P',P_pump(1),'H',h_pump(1),fluid);

for j=2:n_subd
    h_s = CoolProp.PropsSI('H','P',P_pump(j),'S',s_pump(j-1),fluid);
    h_pump(j) = pump('h_in',h_pump(j-1),'h_s',h_s,'eta_p',eta_p);
    T_pump(j) = CoolProp.PropsSI('T','P',P_pump(j),'H',h_pump(j),fluid);
    s_pump(j) = CoolProp.PropsSI('S','P',P_pump(j),'H',h_pump(j),fluid);
end

h_s = CoolProp.PropsSI('H','P',P2,'S',s_pump(1),fluid);
h2 = pump('h_in',h1,'h_s',h_s,'eta_p',eta_p);
T2 = CoolProp.PropsSI('T','P',P2,'H',h2,fluid);
s2 = CoolProp.PropsSI('S','P',P2,'H',h2,fluid);


figure('Color',[1 1 1])
plot(s_pump./1000,T_pump)
hold on
plot(s_pump(1)/1000,T_pump(1),'ro')
plot(s2/1000,T2,'ro')