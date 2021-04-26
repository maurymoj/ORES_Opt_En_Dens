function C_pump = pump_cost(W_p_dot,P,Type,Mat)
% PUMP_COST Estimates pump cost as a function of power, pressure, type of
% pump and material
%   pump_cost()
% Inputs:
%   W_p_dot - Pump power [kW]
%   P       - Pressure   [Pa]

% disp('P_atm = 101.325 kPa and cost index for year 2018.')

Cost_index_ratio = (cost_index(2019)/cost_index(2001));  % Reference 2001 because data from Turton was gathered during the summer of this year

P_bar = P./100000;
P_atm = 101.325/100;
P_barg = P_bar - P_atm;

switch Type
    case 'Reciprocating'
        if (W_p_dot < 0.1) || (W_p_dot > 200)
            error('Power outside valid region');
        elseif P_barg > 100
            error('Pressure exceeds valid region');
        end
        
        K_p = [ 3.8696 0.3161  0.1220];
        
        if P_barg < 10
            C_p = [0 0 0];
        else
            C_p = [-0.245382 0.259016 -0.01363];
        end
        
        switch Mat
            case 'Cast Iron'
                F_m_pump = 1.0;
            case 'Carbon Steel'
                F_m_pump = 1.4;
            case 'Cu alloy'
                F_m_pump = 1.3;
            case 'Stainless Steel'
                F_m_pump = 2.4;
            case 'Ni alloy'
                F_m_pump = 4.0;
            case 'Ti'
                F_m_pump = 6.4;
            otherwise
                F_m_pump = 1.4;
        end
        
    case 'Pos_disp'
        if (W_p_dot < 1) || (W_p_dot > 100)
            error('Power outside valid region');
        elseif P_barg > 100
            error('Pressure exceeds valid region');
        end
        
        K_p = [ 3.4771 0.1350  0.1438];
        
        if P_barg < 10
            C_p = [0 0 0];
        else
            C_p = [-0.245382 0.259016 -0.01363];
        end
        
        switch Mat
            case 'Cast Iron'
                F_m_pump = 1.0;
            case 'Carbon Steel'
                F_m_pump = 1.4;
            case 'Cu alloy'
                F_m_pump = 1.2;
            case 'Stainless Steel'
                F_m_pump = 2.6;
            case 'Ni alloy'
                F_m_pump = 4.8;
            case 'Ti'
                F_m_pump = 10.6;
            otherwise
                F_m_pump = 1.4;
        end
    case 'Centrifugal'
        if (W_p_dot < 1) || (W_p_dot > 300)
            error('Power outside valid region');
        elseif P_barg > 100
            error('Pressure exceeds valid region');
        end
        
        K_p = [ 3.3892 0.0536  0.1538];
        
        if P_barg < 10
            C_p = [0 0 0];
        else
            C_p = [-0.3935 0.3957 -0.00226];
        end
        
        switch Mat
            case 'Cast Iron'
                F_m_pump = 1.0;
            case 'Carbon Steel'
                F_m_pump = 1.6;
            case 'Stainless Steel'
                F_m_pump = 2.2;
            case 'Ni alloy'
                F_m_pump = 4.4;
            otherwise
                F_m_pump = 1.6;
        end
    otherwise
        if (W_p_dot < 0.1) || (W_p_dot > 200)
            error('Power outside valid region');
        elseif P_barg > 100
            error('Pressure exceeds valid region');
        end
        
        K_p = [ 3.3892 0.0536  0.1538];
        
        if P_barg < 10
            C_p = [0 0 0];
        else
            C_p = [-0.3935 0.3957 -0.00226];
        end
        
        switch Mat
            case 'Cast Iron'
                F_m_pump = 1.0;
            case 'Carbon Steel'
                F_m_pump = 1.6;
            case 'Stainless Steel'
                F_m_pump = 2.2;
            case 'Ni alloy'
                F_m_pump = 4.4;
            otherwise
                F_m_pump = 1.6;
        end
end

B_p = [ 1.8900 1.3500];

C_0_pump = 10^(K_p(1) + K_p(2)*log10(W_p_dot) + K_p(3)*(log10(W_p_dot)).^2);
F_P_pump = 10^(C_p(1) + C_p(2)*log10(P_barg) + C_p(3)*(log10(P_barg)).^2);

C_pump = (Cost_index_ratio * C_0_pump) * (B_p(1) + (B_p(2)*F_m_pump*F_P_pump));

end