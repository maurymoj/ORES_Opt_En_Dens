function [c,ceq] = mycon(x)
    eta_RT = ORES_st('P_HPT',x(3));
    c = 0.5 - eta_RT;
    ceq = [];
end