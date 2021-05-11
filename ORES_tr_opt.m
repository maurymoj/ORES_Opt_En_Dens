function [Obj] = ORES_tr_opt(x)

    K_V_H = x(1);
    K_V_L = x(2);
    P_H = x(3);
    
    [eta_RT,rho_E] = ORES_tr('K_V_HPT',K_V_H,'K_V_LPT',K_V_L,'P_HPT',P_H);

%     Obj = -rho_E;
%     Obj = -eta_RT;
    Obj = -rho_E;
end