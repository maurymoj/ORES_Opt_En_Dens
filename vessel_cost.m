function Cost_vessel = vessel_cost(Method,P,Volume,Type,Mat)
% VESSEL_COST Estimates process vessel cost as a function of pressure, type of
% vessel and material
%   vessel_cost()
% Inputs:
%   Method  - 'Bare Module' or 'Matheus'
%   P       - Pressure   [Pa]
%   V       - Volume     [m3]
switch Method
    case 'Bare Module'
        warning('P_atm = 101.325 kPa and cost index for year 2019.')

        P_atm = 101325/100000; % atmospheric pressur in bar
        P_barg = P./100000 - P_atm;

        Cost_index_ratio = (cost_index(2019)/cost_index(2001));  % Reference 2001 because data from Turton was gathered during the summer of this year

        switch Type
            % Turton et al suggests Horizontal vessels for volumes 3.8 < V < 38 m^3
            % and vertical tanks for V < 3.8 and V > 38 m^3
            case 'Horizontal'
                if Volume < 0.1
                    error('Volume too small for cost correlation (smaller than 0.1 m^3).')
                elseif Volume > 628
                    error('Volume too high for cost correlation (higher than 628 m^3).')
                end

                K_v = [3.5565 0.3776 0.0905];
                B_v = [1.49 1.52];

            case 'Vertical'
                if Volume < 0.3
                    error('Volume too small for cost correlation (smaller than 0.3 m^3).')
                elseif Volume > 520
                    error('Volume too high for cost correlation (higher than 520 m^3).')
                end
                K_v = [3.4974 0.4485 0.1074];
                B_v = [2.25 1.82];

            otherwise
                if Volume < 0.3
                    error('Volume too small for cost correlation (smaller than 0.3 m^3).')
                elseif Volume > 520
                    error('Volume too high for cost correlation (higher than 520 m^3).')
                end
                K_v = [3.4974 0.4485 0.1074];
                B_v = [2.25 1.82];
        end

        switch Mat
            case 'Carbon Steel'
                F_m_vessel = 1.0;
                S = 1290; % bar - Turton et al. 2018
        %     case 'Stainless Steel clad'
        %         F_m_vessel = 1.8;
        %         S = 1290;
            case 'Stainless Steel'
                F_m_vessel = 3.1;
                S = 779;
        %     case 'Ni alloy clad'
        %         F_m_vessel = 3.6;
        %         S = 689;
            case 'Ni alloy'
                F_m_vessel = 7.0;
                S = 689;
        %     case 'Ti clad'
        %         F_m_vessel = 4.6;
        %         S = 1041;
            case 'Ti'
                F_m_vessel = 9.4;
                S = 393;
            otherwise
                F_m_vessel = 3.1;
                S = 1290;
        end

        C_0_vessel = 10.^(K_v(1) + K_v(2)*log10(Volume) + K_v(3)*(log10(Volume)).^2);

        E = 0.9; % Weld efficiency - 0.9 corresponds to Single-Welded Butt Joints with Backing Strip

        % Corrosion allowance 
        % CA = 0.0089;  % Known corrosive conditions
        CA = 0.0038;  % noncorrosive streams
        % CA = 0.0015;    % Steam drums and air receivers

        D = (3*Volume/(4*pi))^(1/3); %  Pressure Vessel Design Manual - 2013 - pgs 96 to 98
        % t_min = 0.004763; % Minimum thickness based on ASME
        if D <= 1.07 % minimum thickness based on internal diameter
            t_min = 0.0064; 
        elseif D > 1.07 && D <= 1.52
            t_min = 0.0081;
        elseif D > 1.52
            t_min = 0.0117;
        end

        F_P_vessel = ( P_barg.*D./(2*S*E - 1.2.*P_barg) + CA )./t_min;
        if P_barg < -0.5
            F_P_vessel = 1.25;
        elseif F_P_vessel < 1
            F_P_vessel = 1;
        end

        Cost_vessel = (Cost_index_ratio * C_0_vessel) * (B_v(1) + (B_v(2)*F_m_vessel*F_P_vessel));
    case 'Matheus'
        Storage_Cost_per_m3 = 750; % Pressure tank cost [$/m3]
        Cost_vessel = Storage_Cost_per_m3*Volume;
end

end