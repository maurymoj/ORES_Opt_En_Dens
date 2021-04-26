function C_turb = turbine_cost(Method,W_t_dot,Type,Mat)
% TURBINE_COST Estimates turbine cost as a function of power, type of
% turbine and material
%   turbine_cost()
% Inputs:
%   W_t_dot - Turbine power [kW]
% Methods:
%   Bare Module - from Turton et al. 2018
%   six-tenth rule

switch Method
    case 'Bare Module'
%         warning('P_atm = 101.325 kPa and cost index for year 2018.')

        Cost_index_ratio = (cost_index(2018)/cost_index(2001));  % Reference 2001 because data from Turton was gathered during the summer of this year

%         P_bar = P./100000;
%         P_atm = 101.325/100;
%         P_barg = P_bar - P_atm;

        switch Type
            case 'Axial Gas Turbine'
                if (W_t_dot < 100) || (W_t_dot > 4000)
                    error('Power outside valid region');
                end

                K_t = [ 2.7051 1.4398  -0.1776];

                switch Mat
                    case 'Carbon Steel'
                        F_m_turb = 3.5;
                    case 'Stainless Steel'
                        F_m_turb = 6.2;
                    case 'Ni alloy'
                        F_m_turb = 11.6;           
                    otherwise
                        F_m_turb = 3.5;
                end

            case 'Radial gas-liquid expander'
                if (W_t_dot < 100) || (W_t_dot > 1500)
                    error('Power outside valid region');
                end

                K_t = [ 2.2476 1.4965 -0.1618];

                switch Mat
                    case 'Carbon Steel'
                        F_m_turb = 3.5;
                    case 'Stainless Steel'
                        F_m_turb = 6.2;
                    case 'Ni alloy'
                        F_m_turb = 11.6;           
                    otherwise
                        F_m_turb = 3.5;
                end

            otherwise
                if (W_t_dot < 0.1) || (W_t_dot > 200)
                    error('Power outside valid region');
                end

                K_t = [ 2.7051 1.4398  -0.1776];

                switch Mat
                    case 'Carbon Steel'
                        F_m_turb = 3.5;
                    case 'Stainless Steel'
                        F_m_turb = 6.2;
                    case 'Ni alloy'
                        F_m_turb = 11.6;           
                    otherwise
                        F_m_turb = 3.5;
                end
        end

        C_0_turb = 10^(K_t(1) + K_t(2)*log10(W_t_dot) + K_t(3)*(log10(W_t_dot)).^2);
        C_turb = (Cost_index_ratio * C_0_turb) * F_m_turb;
    case 'six-tenth'
        Cost_index_ratio = (cost_index(2018)/cost_index(2009));  % Reference 2009 because data from ESMAP Technical Paper 122/09
        C_turb = Cost_index_ratio*(130*W_t_dot)*(W_t_dot/325000).^0.6;
end

end