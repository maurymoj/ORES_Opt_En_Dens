function [Obj] = ORES_st_opt(x)
    P_H = x(1);
    Delta_T = x(2);
    fluid = 'R141b';
    
    T = CoolProp.PropsSI('T','P',P_H,'Q',1,fluid)+Delta_T;
    
    [eta_RT,eta_I,eta_II,q_in,w_net]=ORES_st('P_HPT',P_H,'T_2',T,'fluid',fluid);
    
    Obj = -eta_RT;
end