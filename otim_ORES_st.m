function [eta_RT] = otim_ORES_st(x)
% otim_ORES_st
%   Maximização da eta_RT em função de P_HPT e DT_SH para operação em
%   regime permanente
P = x(1)*1e3;
DT = x(2);

[eta_I,eta_II,w_net,eta_RT] = ORES_st('fluid','R134a','P_HPT',P,'DT_SH',DT);
eta_RT = -eta_RT;
end