V = 1;
V_HPT = V;
V_LPT = V;

% P_HPT = 2000000;
P_HPT=1000000:250000:3000000;
P_LPT = 600000:100000:900000;

Q_LPT_f = zeros(length(P_HPT),1);

for j=1:length(P_LPT)
    for i=1:length(P_HPT)
        rho_HPT_0=CoolProp.PropsSI('D','P',P_HPT(i),'Q',0,'R134a');
        rho_HPT_1=CoolProp.PropsSI('D','P',P_HPT(i),'Q',1,'R134a');
        rho_LPT_1=CoolProp.PropsSI('D','P',P_LPT(j),'Q',1,'R134a');

        m_HPT_i = rho_HPT_0*V_HPT;
        m_HPT_f = rho_HPT_1*V_HPT;
        Dm_HPT = m_HPT_i - m_HPT_f;

        m_LPT_i = rho_LPT_1*V;
        m_LPT_f = m_LPT_i + Dm_HPT;
        rho_LPT_f = m_LPT_f/V;
        Q_LPT_f(i,j) = CoolProp.PropsSI('Q','P',P_LPT(j),'D',rho_LPT_f,'R134a');
    end
end

figure('color',[1 1 1])
hold all
grid on

for j=1:length(P_LPT)
    plot(P_HPT./1000,Q_LPT_f(:,j))
end

xlabel('P_{HPT} [kPa]')
ylabel('x_{LPT,f}')
legend(num2str(P_LPT'))