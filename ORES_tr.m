function [eta_RT,rho_E_L,t_generation,erro,CAPEX]= ORES_tr(varargin)
% ORES_TR Transient operation model
%   ORES_tr()

% sets of inputs and outputs:
%=========================================================================%    
% Como configurar referência para a propriedades termodinâmicas:          %
% CoolProp.set_reference_stateS('R134a','ASHRAE')                         %
%=========================================================================%

%-------------------------- Simulation outputs ---------------------------%
    if any(strcmp('Plot',varargin))
        plotBool = varargin{find(strcmp('Plot',varargin))+1};
        if ~isnumeric(plotBool)
            error('Plot must be a number or boolean, not a %s.',class(plotBool))
        end
    else
        plotBool = false;
    end
    
%====================== Assignment of input values =======================%

%--------------------------- Default parameters --------------------------%
% Ambient conditions
T_amb = 298.15;
P_amb = 101325;

% System parameters and requirements
fluid = 'R141b';
% W_g = 2000;     % Generated power
W_g = 1000000;  % 1 MW
W_c = 0.2*W_g;% Consumed power

% W_t = 2000;     % Turbine power
% W_p = 0.2*W_t;  % Pump power

% Simulation parameters
dt = 25;     % [s] 50 s also provides good results (about 1% higher than with 25 s
Delta_t = 3600; % Total energy generation time [s]

% Assumptions
eta_t = 0.80;
eta_p = 0.75;

P_HPT_default = 3200000; % 3.2 MPa
x_HPT_default = 0.02;

% P_LPT_default = 1000000;
T_LPT_default = T_amb;
x_LPT_default = 0.98;
DT_SH = 0.01;

x_HPT_f_disch = 0.98;
x_LPT_f_disch = 0.02;

x_HPT_f_ch = 0.02;
x_LPT_f_ch = 0.98;

%
K_V_H = 2; % Ratio of volume of the HPT over the minimum required volume
K_V_L = 1; % Ratio of V_LPT over V_HPT

Heat_source = 'HPs';
COP_HP1_0 = 0;
COP_HP2_0 = 0;

erro = 0;

%------- Basic parameters: Fluid definition and ambient conditions -------%
    if any(strcmp('fluid',varargin))
        fluid = varargin{find(strcmp('fluid',varargin))+1};
        if ~ischar(fluid)
            error('fluid must be a string, not a %s.',class(fluid))
        end
    end
%     Fluid_Cost_per_kg = [4.49 7.01 4.41 2.00 11.50 5.00]; % Cost based on Alibaba search on 16/01/2020 - Matheus
    switch fluid
        case 'R152a'
            rho_lim = [360,375];
            P_lim_inf = 1800000;
            %Fluid_Cost_per_kg = 4.49;
        case 'R134a'
            rho_lim = [443,514];    % INSTÁVEL
            P_lim_inf = 1700000;
            %Fluid_Cost_per_kg = 7.01;
        case 'R142b'
            rho_lim = [0 0];
            P_lim_inf = 4.41;
            %Fluid_Cost_per_kg = 4.41;
        case 'R365mfc'
            rho_lim = [471,491];
            P_lim_inf = 2.00;
            %Fluid_Cost_per_kg = 2.00;
        case 'R236ea'
            rho_lim = [0 0];
            P_lim_inf = 800000;
            %Fluid_Cost_per_kg = 6.00; % Valor arbitrado - não foi encontrado valor no Alibaba
        case 'R141b'
            rho_lim = [0 0];
            P_lim_inf = 700000;
            %Fluid_Cost_per_kg = 5.00;
        case 'CO2'
            rho_lim = [0 0];
            P_lim_inf = 640000;
        otherwise
            rho_lim = [0,0];
            P_lim_inf = 0;
    end
    
    if any(strcmp('T_amb',varargin))
        T_amb = varargin{find(strcmp('T_amb',varargin))+1};
        if ~isnumeric(T_amb)
            error('T_amb must be a number, not a %s.',class(T_amb))
        end
    end
    if any(strcmp('P_amb',varargin))
        P_amb = varargin{find(strcmp('P_amb',varargin))+1};
        if ~isnumeric(P_amb)
            error('P_amb must be a number, not a %s.',class(P_amb))
        end
    end
    % Dead state properties
    P_o = P_amb;
    T_o = T_amb;
    u_o = CoolProp.PropsSI('U','P',P_o,'T',T_o,fluid);
%     h_o = CoolProp.PropsSI('H','P',P_o,'T',T_o,fluid);
    v_o = 1/CoolProp.PropsSI('D','P',P_o,'T',T_o,fluid);
    s_o = CoolProp.PropsSI('S','P',P_o,'T',T_o,fluid);
    rho_o = CoolProp.PropsSI('D','P',P_o,'T',T_o,fluid);
    
%---------------------------- Power generation ---------------------------%
    if any(strcmp('W_g',varargin))
        W_g = varargin{find(strcmp('W_g',varargin))+1};
        if ~isnumeric(W_g)
            error('W_g must be a number, not a %s.',class(W_g))
        end
    end  
%---------------------------- Power consumption --------------------------%
    if any(strcmp('W_c',varargin))
        W_c = varargin{find(strcmp('W_c',varargin))+1};
        if ~isnumeric(W_c)
            error('W_c must be a number, not a %s.',class(W_c))
        end
    end

%------------------------------- Heat Source -----------------------------%
    if any(strcmp('Heat_Source',varargin))
        Heat_source = varargin{find(strcmp('Heat_Source',varargin))+1};
        if isnumeric(Heat_source)
            error('Heat_source must be a string, not a %s.',class(Heat_source))
        end
        switch Heat_source
            case 'HPs'      % Heat pumps
                if any(strcmp('COP_HP',varargin))
                    COP_HP1_0 = varargin{find(strcmp('COP_HP',varargin))+1};
                    COP_HP2_0 = varargin{find(strcmp('COP_HP',varargin))+1};
                    if ~isnumeric(COP_HP1_0)
                        error('COP_HP must be a number, not a %s.',class(COP_HP1_0))
                    end
                elseif any(strcmp('COP_HP1',varargin))
                    COP_HP1_0 = varargin{find(strcmp('COP_HP1',varargin))+1};
                    if ~isnumeric(COP_HP1_0)
                        error('COP_HP1 must be a number, not a %s.',class(COP_HP1_0))
                    end
                    if any(strcmp('COP_HP2',varargin))
                        COP_HP2_0 = varargin{find(strcmp('COP_HP2',varargin))+1};
                        if ~isnumeric(COP_HP2_0)
                            error('COP_HP2 must be a number, not a %s.',class(COP_HP2_0))
                        end
                    end
                elseif any(strcmp('COP_HP2',varargin))
                    COP_HP2_0 = varargin{find(strcmp('COP_HP2',varargin))+1};
                    if ~isnumeric(COP_HP2_0)
                        error('COP_HP2 must be a number, not a %s.',class(COP_HP2_0))
                    end                  
                end
            case 'Free'     % Heat is provided for free
            % Default - Reversible heat pumps 
        end

    end
%======================= READING INPUT PARAMETERS ========================%

%------------------------ SIMULATION PARAMETERS --------------------------%
%--- Time step, Total discharging time and number of steps during --------%
%--- discharge -----------------------------------------------------------%

    if any(strcmp('dt',varargin)) 
        dt = varargin{find(strcmp('dt',varargin))+1};
        if ~isnumeric(dt)
            error('dt must be a number, not a %s.',class(dt))
        end
        
        if any(strcmp('Delta_t',varargin)) % If dt and Delta_t are given.
            Delta_t = varargin{find(strcmp('Delta_t',varargin))+1};
            if ~isnumeric(Delta_t)
                error('Delta_t must be a number, not a %s.',class(Delta_t))
            end
            N_Steps = floor(Delta_t/dt);
            if N_Steps<1
                error('dt must be smaller than Delta_t.')
            end
        elseif any(strcmp('N_Steps',varargin)) % If dt and N_Steps are given.
            N_Steps = varargin{find(strcmp('N_Steps',varargin))+1};
            if ~isnumeric(N_Steps)
                error('N_Steps must be a number, not a %s.',class(N_Steps))
            elseif N_Steps<1
                error('N_Steps must be at least equal to 1.')
            end
            N_Steps = floor(N_Steps);
            Delta_t = dt*N_Steps;
        else                                % If only dt is given.
            N_Steps = floor(Delta_t/dt);
            if N_Steps<1
                error('dt must be smaller than Delta_t.')
            end
        end
    elseif any(strcmp('Delta_t',varargin))
        Delta_t = varargin{find(strcmp('Delta_t',varargin))+1};
        if ~isnumeric(Delta_t)
            error('Delta_t must be a number, not a %s.',class(Delta_t))
        end
        if any(strcmp('N_Steps',varargin)) % If Delta_t and N_Steps are given.
            N_Steps = varargin{find(strcmp('N_Steps',varargin))+1};
            if ~isnumeric(N_Steps)
                error('N_Steps must be a number, not a %s.',class(N_Steps))
            elseif N_Steps<1
                error('N_Steps must be at least equal to 1.')
            end
            dt = Delta_t/N_Steps;
        else                               % If only Delta_t is given.
            N_Steps = floor(Delta_t/dt);
            if N_Steps<1
                error('Delta_t must be larger than dt.')
            end
        end 
    elseif any(strcmp('N_Steps',varargin)) % If only N_Steps is given.
        N_Steps = varargin{find(strcmp('N_Steps',varargin))+1};
        if ~isnumeric(N_Steps)
            error('N_Steps must be a number, not a %s.',class(N_Steps))
        elseif N_Steps<1
                error('N_Steps must be at least equal to 1.')
        end
        N_Steps = floor(N_Steps);
        Delta_t = dt*N_Steps;
    else                                   % DEFAULT
        N_Steps = floor(Delta_t/dt);
    end
    
%------------ Initialization of discharging phase parameters -------------%
    t_disch = zeros(N_Steps+1,1);
    P_LPT = zeros(N_Steps+1,1);
    x_LPT = zeros(N_Steps+1,1);
    T_LPT = zeros(N_Steps+1,1);
    rho_LPT = zeros(N_Steps+1,1);
    u_LPT = zeros(N_Steps+1,1);
    h_LPT = zeros(N_Steps+1,1);
    s_LPT = zeros(N_Steps+1,1);
    ex_LPT = zeros(N_Steps+1,1);
    m_LPT = zeros(N_Steps+1,1);
    U_LPT = zeros(N_Steps+1,1);
    EX_LPT = zeros(N_Steps+1,1);
    
    P_HPT = zeros(N_Steps+1,1);
    x_HPT = zeros(N_Steps+1,1);
    T_HPT = zeros(N_Steps+1,1);
    rho_HPT = zeros(N_Steps+1,1);
%     ph_HPT = zeros(N_Steps+1,1);    
    u_HPT = zeros(N_Steps+1,1);
    h_HPT = zeros(N_Steps+1,1);
    s_HPT = zeros(N_Steps+1,1);
    ex_HPT = zeros(N_Steps+1,1);
    m_HPT = zeros(N_Steps+1,1);
    U_HPT = zeros(N_Steps+1,1);
    EX_HPT = zeros(N_Steps+1,1);   
    
    P_1 = zeros(N_Steps+1,1);
    T_1 = zeros(N_Steps+1,1);
    s_1 = zeros(N_Steps+1,1);
    h_1 = zeros(N_Steps+1,1);
    T_sat_Ev = zeros(N_Steps+1,1);
    
    P_2 = zeros(N_Steps+1,1);
    T_2 = zeros(N_Steps+1,1);
    s_2 = zeros(N_Steps+1,1);
    h_2 = zeros(N_Steps+1,1);
    
    P_3 = zeros(N_Steps+1,1);
    T_3 = zeros(N_Steps+1,1);
    s_3 = zeros(N_Steps+1,1);
    h_t = zeros(N_Steps+1,1);
    h_3 = zeros(N_Steps+1,1);
    
    P_4 = zeros(N_Steps+1,1);
    T_4 = zeros(N_Steps+1,1);
    s_4 = zeros(N_Steps+1,1);
    h_4 = zeros(N_Steps+1,1);
    x_4 = zeros(N_Steps+1,1);
    
    q_ev = zeros(N_Steps+1,1);
    w_t = zeros(N_Steps+1,1);
    w_g = zeros(N_Steps+1,1);
            
    m_disch = zeros(N_Steps+1,1);
%     dE_g = zeros(N_Steps+1,1);
    dE_t = zeros(N_Steps+1,1);
    dQ_ev = zeros(N_Steps+1,1);
    
    dE_net = zeros(N_Steps+1,1);

    switch Heat_source
        case 'HPs'
            h_g_ev = zeros(N_Steps+1,1);
            q_SH = zeros(N_Steps+1,1);
            q_lg = zeros(N_Steps+1,1);
            T_H1 = zeros(N_Steps+1,1);
            COP_carnot_HP1 = zeros(N_Steps+1,1);
            COP_HP1 = zeros(N_Steps+1,1);
            w_c_HP1 = zeros(N_Steps+1,1);
            q_L_HP1 = zeros(N_Steps+1,1);
            dQ_L_HP1 = zeros(N_Steps+1,1);
        otherwise
            h_g_ev = zeros(N_Steps+1,1);
            q_SH = zeros(N_Steps+1,1);
            q_lg = zeros(N_Steps+1,1);
            T_H1 = zeros(N_Steps+1,1);
            COP_carnot_HP1 = zeros(N_Steps+1,1);
            COP_HP1 = zeros(N_Steps+1,1);
            w_c_HP1 = zeros(N_Steps+1,1);
            q_L_HP1 = zeros(N_Steps+1,1);
            dQ_L_HP1 = zeros(N_Steps+1,1);
%            dE_HP1 = zeros(N_Steps+1,1);
    end

%-------------------------- Initial state at HPT -------------------------%
    if any(strcmp('P_HPT',varargin))    % INPUT P_HPT
        if ~isnumeric(varargin{find(strcmp('P_HPT',varargin))+1})
            error('P_HPT must be a number, not a %s.',class(varargin{find(strcmp('P_HPT',varargin))+1}))
        end
        P_HPT(1) = varargin{find(strcmp('P_HPT',varargin))+1};
        if any(strcmp('x_HPT',varargin))    % INPUT P_HPT + x_HPT
            if ~isnumeric(varargin{find(strcmp('x_HPT',varargin))+1})
                error('x_HPT must be a number, not a %s.',class(varargin{find(strcmp('x_HPT',varargin))+1}))
            elseif any(strcmp('T_HPT',varargin))
                warning('P_HPT and x_HPT are already determined, given T_HPT value will be ignored.');
            end
            x_HPT(1) = varargin{find(strcmp('x_HPT',varargin))+1};
            T_HPT(1) = CoolProp.PropsSI('T','P',P_HPT(1),'Q',x_HPT(1),fluid);
            rho_HPT(1) = CoolProp.PropsSI('D','P',P_HPT(1),'Q',x_HPT(1),fluid);
            u_HPT(1) = CoolProp.PropsSI('U','P',P_HPT(1),'Q',x_HPT(1),fluid);
            h_HPT(1) = CoolProp.PropsSI('H','P',P_HPT(1),'Q',x_HPT(1),fluid);
            s_HPT(1) = CoolProp.PropsSI('S','P',P_HPT(1),'Q',x_HPT(1),fluid);
            ph_HPT(1) = CoolProp.PropsSI('Phase','P',P_HPT(1),'Q',x_HPT(1),fluid);
        elseif any(strcmp('T_HPT',varargin))    % INPUT P_HPT + T_HPT
            if ~isnumeric(varargin{find(strcmp('T_HPT',varargin))+1})
                error('T_HPT must be a number, not a %s.',class(varargin{find(strcmp('T_HPT',varargin))+1}))
            end
            T_HPT(1) = varargin{find(strcmp('T_HPT',varargin))+1};
            %x_HPT(1) = CoolProp.PropsSI('Q','P',P_HPT,'T',T_HPT,fluid);
            rho_HPT(1) = CoolProp.PropsSI('D','P',P_HPT(1),'T',T_HPT(1),fluid);
            u_HPT(1) = CoolProp.PropsSI('U','P',P_HPT(1),'T',T_HPT(1),fluid);
            h_HPT(1) = CoolProp.PropsSI('H','P',P_HPT(1),'T',T_HPT(1),fluid);
            s_HPT(1) = CoolProp.PropsSI('S','P',P_HPT(1),'T',T_HPT(1),fluid);
            ph_HPT(1) = CoolProp.PropsSI('Phase','P',P_HPT(1),'T',T_HPT(1),fluid);
        else                % INPUT P_HPT + DEFAULT
            x_HPT(1) = x_HPT_default;
            T_HPT(1) = CoolProp.PropsSI('T','P',P_HPT(1),'Q',x_HPT(1),fluid);
            rho_HPT(1) = CoolProp.PropsSI('D','P',P_HPT(1),'Q',x_HPT(1),fluid);
            u_HPT(1) = CoolProp.PropsSI('U','P',P_HPT(1),'Q',x_HPT(1),fluid);
            h_HPT(1) = CoolProp.PropsSI('H','P',P_HPT(1),'Q',x_HPT(1),fluid);
            s_HPT(1) = CoolProp.PropsSI('S','P',P_HPT(1),'Q',x_HPT(1),fluid);
            ph_HPT(1) = CoolProp.PropsSI('Phase','P',P_HPT(1),'Q',x_HPT(1),fluid);
        end
    else
        if any(strcmp('x_HPT',varargin))        % INPUT x_HPT
            if ~isnumeric(varargin{find(strcmp('x_HPT',varargin))+1})
                error('x_HPT must be a number, not a %s.',class(varargin{find(strcmp('x_HPT',varargin))+1}))
            end
            x_HPT(1) = varargin{find(strcmp('x_HPT',varargin))+1};
            if any(strcmp('T_HPT',varargin))    % INPUT x_HPT + T_HPT
                if ~isnumeric(varargin{find(strcmp('T_HPT',varargin))+1})
                    error('T_HPT must be a number, not a %s.',class(varargin{find(strcmp('T_HPT',varargin))+1}))
                end
                T_HPT(1) = varargin{find(strcmp('T_HPT',varargin))+1};
                P_HPT(1) = CoolProp.PropsSI('P','Q',x_HPT(1),'T',T_HPT(1),fluid);
                rho_HPT(1) = CoolProp.PropsSI('D','Q',x_HPT(1),'T',T_HPT(1),fluid);
                u_HPT(1) = CoolProp.PropsSI('U','Q',x_HPT(1),'T',T_HPT(1),fluid);
                h_HPT(1) = CoolProp.PropsSI('H','Q',x_HPT(1),'T',T_HPT(1),fluid);
                s_HPT(1) = CoolProp.PropsSI('S','Q',x_HPT(1),'T',T_HPT(1),fluid);
                ph_HPT(1) = CoolProp.PropsSI('Phase','Q',x_HPT(1),'T',T_HPT(1),fluid);
            else                            % INPUT x_HPT + DEFAULT
                P_HPT(1) = P_HPT_default;
                T_HPT(1) = CoolProp.PropsSI('T','Q',x_HPT(1),'P',P_HPT(1),fluid);
                rho_HPT(1) = CoolProp.PropsSI('D','Q',x_HPT(1),'P',P_HPT(1),fluid);
                u_HPT(1) = CoolProp.PropsSI('U','Q',x_HPT(1),'P',P_HPT(1),fluid);
                h_HPT(1) = CoolProp.PropsSI('H','Q',x_HPT(1),'P',P_HPT(1),fluid);
                s_HPT(1) = CoolProp.PropsSI('S','Q',x_HPT(1),'P',P_HPT(1),fluid);
                ph_HPT(1) = CoolProp.PropsSI('Phase','Q',x_HPT(1),'P',P_HPT(1),fluid);
            end
        elseif any(strcmp('T_HPT',varargin))  % INPUT T_HPT + DEFAULT
            if ~isnumeric(varargin{find(strcmp('T_HPT',varargin))+1})
                error('T_HPT must be a number, not a %s.',class(varargin{find(strcmp('T_HPT',varargin))+1}))
            end
            T_HPT(1) = varargin{find(strcmp('T_HPT',varargin))+1};
            x_HPT(1) = x_HPT_default;
            P_HPT(1) = CoolProp.PropsSI('P','Q',x_HPT(1),'T',T_HPT(1),fluid);
            rho_HPT(1) = CoolProp.PropsSI('D','Q',x_HPT(1),'T',T_HPT(1),fluid);
            u_HPT(1) = CoolProp.PropsSI('U','Q',x_HPT(1),'T',T_HPT(1),fluid);
            h_HPT(1) = CoolProp.PropsSI('H','Q',x_HPT(1),'T',T_HPT(1),fluid);
            s_HPT(1) = CoolProp.PropsSI('S','Q',x_HPT(1),'T',T_HPT(1),fluid);
            ph_HPT(1) = CoolProp.PropsSI('Phase','Q',x_HPT(1),'T',T_HPT(1),fluid);
        else                                % If no variable for 4 is available
%             warning('No input for HPT, default values for x=0 and P=3.500 kPa assumed.')
            P_HPT(1) = P_HPT_default;
            x_HPT(1) = x_HPT_default;
            T_HPT(1) = CoolProp.PropsSI('T','Q',x_HPT(1),'P',P_HPT(1),fluid);
            rho_HPT(1) = CoolProp.PropsSI('D','Q',x_HPT(1),'P',P_HPT(1),fluid);
            u_HPT(1) = CoolProp.PropsSI('U','Q',x_HPT(1),'P',P_HPT(1),fluid);
            h_HPT(1) = CoolProp.PropsSI('H','Q',x_HPT(1),'P',P_HPT(1),fluid);
            s_HPT(1) = CoolProp.PropsSI('S','Q',x_HPT(1),'P',P_HPT(1),fluid);
            ph_HPT(1) = CoolProp.PropsSI('Phase','Q',x_HPT(1),'P',P_HPT(1),fluid);
        end
    end
    ex_HPT(1) = (u_HPT(1)-u_o) - T_o*(s_HPT(1)-s_o) + P_o*(1./rho_HPT(1)-1/rho_o);
    
%------------ PRESSURE DROP IN THE EVAPORATOR/SUPERHEATER ----------------%
    if any(strcmp('DP_Ev',varargin))
        DeltaP_Ev = varargin{find(strcmp('DP_Ev',varargin))+1};
        if ~isnumeric(varargin{find(strcmp('DP_Ev',varargin))+1})
            error('DP_Ev must be a number, not a %s.',class(varargin{find(strcmp('DP_Ev',varargin))+1}))
        end
% ACRESCENTAR FORMAS DE CALCULAR A PERDA DE PRESSÃO EM FUNÇÃO DE OUTROS PARÂMETROS                
    else
        DeltaP_Ev = 0;
    end    

%----------- STATE AT TURBINE INLET AND ISOENTROPIC EFFICIENCY -----------%
    if any(strcmp('T_2',varargin)) && any(strcmp('T_H',varargin))
        T_2 = varargin{find(strcmp('T_2',varargin))+1};
%         T_H = varargin{find(strcmp('T_H',varargin))+1};
        if ~isnumeric(varargin{find(strcmp('T_2',varargin))+1})
            error('T_2 must be a number, not a %s.',class(varargin{find(strcmp('T_2',varargin))+1}))
        elseif ~isnumeric(varargin{find(strcmp('T_H',varargin))+1})
            error('T_H must be a number, not a %s.',class(varargin{find(strcmp('T_H',varargin))+1}))
        end        
    elseif any(strcmp('T_2',varargin)) 
        if ~isnumeric(varargin{find(strcmp('T_2',varargin))+1})
            error('T_2 must be a number, not a %s.',class(varargin{find(strcmp('T_2',varargin))+1}))
        end
        T_2(1) = varargin{find(strcmp('T_2',varargin))+1};
        T_sat_Ev(1) = CoolProp.PropsSI('T','P',P_HPT(1),'Q',0,fluid);
%         T_H = T_2(1) + 20;
%     elseif any(strcmp('T_H',varargin))
%         if ~isnumeric(varargin{find(strcmp('T_H',varargin))+1})
%             error('T_H must be a number, not a %s.',class(varargin{find(strcmp('T_H',varargin))+1}))
%         end
%         T_H = varargin{find(strcmp('T_H',varargin))+1};
%         T_2(1) = T_H-20;
    elseif any(strcmp('DT_SH',varargin))
        if ~isnumeric(varargin{find(strcmp('DT_SH',varargin))+1})
            error('DT_SH must be a number, not a %s.',class(varargin{find(strcmp('DT_SH',varargin))+1}))
        end
        DT_SH = varargin{find(strcmp('DT_SH',varargin))+1};
        if DT_SH==0
            warning('DT_SH must be greater than 0. It was adjusted to 0.01 K.')
            DT_SH = 0.01;
        end
        T_sat_Ev(1) = CoolProp.PropsSI('T','P',P_HPT(1),'Q',0,fluid);
        T_2(1) = T_sat_Ev(1) + DT_SH;
%         T_H = T_2(1) + 20;
    else % DEFAULT
        T_sat_Ev(1) = CoolProp.PropsSI('T','P',P_HPT(1),'Q',0,fluid);
        T_2(1) = T_sat_Ev(1) + DT_SH;
%         T_H = T_2(1) + 20;
    end            
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% INCLUIR VARIAÇÃO DA EFICIÊNCIA ISOENTRÓPICA DA TURBINA E DA BOMBA EM
% FUNÇÃO DAS CONDIÇÕES DE OPERAÇÃO
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if any(strcmp('eta_t',varargin))
        if ~isnumeric(varargin{find(strcmp('eta_t',varargin))+1})
            error('eta_t must be a number, not a %s.',class(varargin{find(strcmp('eta_t',varargin))+1}))
        end
        eta_t = varargin{find(strcmp('eta_t',varargin))+1};
    end
    
%------------------- PRESSURE DROP IN THE CONDENSER ---------------------%
    if any(strcmp('DP_Con',varargin))
        DeltaP_Con = varargin{find(strcmp('DP_Con',varargin))+1};
        if ~isnumeric(varargin{find(strcmp('DP_Con',varargin))+1})
            error('DP_Con must be a number, not a %s.',class(varargin{find(strcmp('DP_Con',varargin))+1}))
        end
% ACRESCENTAR FORMAS DE CALCULAR A PERDA DE PRESSÃO EM FUNÇÃO DE OUTROS PARÂMETROS
    else
        DeltaP_Con = 0;
    end
    
%------------------------- Initial state at LPT --------------------------%
    if any(strcmp('P_LPT',varargin))
        if ~isnumeric(varargin{find(strcmp('P_LPT',varargin))+1})
            error('P_LPT must be a number, not a %s.',class(varargin{find(strcmp('P_LPT',varargin))+1}))
        end
        P_LPT(1) = varargin{find(strcmp('P_LPT',varargin))+1};
        if any(strcmp('x_LPT',varargin))      % If P and x are available
            if ~isnumeric(varargin{find(strcmp('x_LPT',varargin))+1})
                error('x_LPT must be a number, not a %s.',class(varargin{find(strcmp('x_LPT',varargin))+1}))
            elseif any(strcmp('T_LPT',varargin))
                warning('P_LPT and x_LPT are already determined, given T_LPT value will be ignored.')
            end
            x_LPT(1) = varargin{find(strcmp('x_LPT',varargin))+1};
            T_LPT(1) = CoolProp.PropsSI('T','P',P_LPT(1),'Q',x_LPT(1),fluid);
            rho_LPT(1) = CoolProp.PropsSI('D','P',P_LPT(1),'Q',x_LPT(1),fluid);
            u_LPT(1) = CoolProp.PropsSI('U','P',P_LPT(1),'Q',x_LPT(1),fluid);
            h_LPT(1) = CoolProp.PropsSI('H','P',P_LPT(1),'Q',x_LPT(1),fluid);
            s_LPT(1) = CoolProp.PropsSI('S','P',P_LPT(1),'Q',x_LPT(1),fluid);
        elseif any(strcmp('T_LPT',varargin))  % If P and x are available
            if ~isnumeric(varargin{find(strcmp('T_LPT',varargin))+1})
                error('T_LPT must be a number, not a %s.',class(varargin{find(strcmp('T_LPT',varargin))+1}))
            end
            T_LPT(1) = varargin{find(strcmp('T_LPT',varargin))+1};
            x_LPT(1) = CoolProp.PropsSI('Q','P',P_LPT(1),'T',T_LPT(1),fluid);
            rho_LPT(1) = CoolProp.PropsSI('D','P',P_LPT(1),'T',T_LPT(1),fluid);
            u_LPT(1) = CoolProp.PropsSI('U','P',P_LPT(1),'T',T_LPT(1),fluid);
            h_LPT(1) = CoolProp.PropsSI('H','P',P_LPT(1),'T',T_LPT(1),fluid);            
            s_LPT(1) = CoolProp.PropsSI('S','P',P_LPT(1),'T',T_LPT(1),fluid);
        else                                % If only P is available
            x_LPT(1) = x_LPT_default;
            T_LPT(1) = CoolProp.PropsSI('T','P',P_LPT(1),'Q',x_LPT(1),fluid);
            rho_LPT(1) = CoolProp.PropsSI('D','P',P_LPT(1),'Q',x_LPT(1),fluid);
            u_LPT(1) = CoolProp.PropsSI('U','P',P_LPT(1),'Q',x_LPT(1),fluid);
            h_LPT(1) = CoolProp.PropsSI('H','P',P_LPT(1),'Q',x_LPT(1),fluid);
            s_LPT(1) = CoolProp.PropsSI('S','P',P_LPT(1),'Q',x_LPT(1),fluid);
        end 
    else
        if any(strcmp('x_LPT',varargin)) 
            if ~isnumeric(varargin{find(strcmp('x_LPT',varargin))+1})
                error('x_LPT must be a number, not a %s.',class(varargin{find(strcmp('x_LPT',varargin))+1}))
            end
            x_LPT(1) = varargin{find(strcmp('x_LPT',varargin))+1};
            if any(strcmp('T_LPT',varargin))  % If x and T are available
                if ~isnumeric(varargin{find(strcmp('T_LPT',varargin))+1})
                    error('T_LPT must be a number, not a %s.',class(varargin{find(strcmp('T_LPT',varargin))+1}))
                end
                T_LPT(1) = varargin{find(strcmp('T_LPT',varargin))+1};
                P_LPT(1) = CoolProp.PropsSI('P','T',T_LPT(1),'Q',x_LPT(1),fluid);
                rho_LPT(1) = CoolProp.PropsSI('D','T',T_LPT(1),'Q',x_LPT(1),fluid);
                u_LPT(1) = CoolProp.PropsSI('U','T',T_LPT(1),'Q',x_LPT(1),fluid);
                h_LPT(1) = CoolProp.PropsSI('H','T',T_LPT(1),'Q',x_LPT(1),fluid);
                s_LPT(1) = CoolProp.PropsSI('S','T',T_LPT(1),'Q',x_LPT(1),fluid);
            else                            % If only x is available
                T_LPT(1) = x_LPT_default;
                P_LPT(1) = CoolProp.PropsSI('P','T',T_LPT(1),'Q',x_LPT(1),fluid);
                rho_LPT(1) = CoolProp.PropsSI('D','T',T_LPT(1),'Q',x_LPT(1),fluid);
                u_LPT(1) = CoolProp.PropsSI('U','T',T_LPT(1),'Q',x_LPT(1),fluid);
                h_LPT(1) = CoolProp.PropsSI('H','T',T_LPT(1),'Q',x_LPT(1),fluid);
                s_LPT(1) = CoolProp.PropsSI('S','T',T_LPT(1),'Q',x_LPT(1),fluid);
            end
        elseif any(strcmp('T_LPT',varargin))  % If only T is available
            if ~isnumeric(varargin{find(strcmp('T_LPT',varargin))+1})
                error('T_LPT must be a number, not a %s.',class(varargin{find(strcmp('T_LPT',varargin))+1}))
            end
            T_LPT(1) = varargin{find(strcmp('T_LPT',varargin))+1};
            x_LPT(1) = x_LPT_default;
            P_LPT(1) = CoolProp.PropsSI('P','T',T_LPT(1),'Q',x_LPT(1),fluid);
            rho_LPT(1) = CoolProp.PropsSI('D','T',T_LPT(1),'Q',x_LPT(1),fluid);
            u_LPT(1) = CoolProp.PropsSI('U','T',T_LPT(1),'Q',x_LPT(1),fluid);
            h_LPT(1) = CoolProp.PropsSI('H','T',T_LPT(1),'Q',x_LPT(1),fluid);
            s_LPT(1) = CoolProp.PropsSI('S','T',T_LPT(1),'Q',x_LPT(1),fluid);
        else                                % If no variable for 4 is available
            T_LPT(1) = T_LPT_default;
            x_LPT(1) = x_LPT_default;
            P_LPT(1) = CoolProp.PropsSI('P','T',T_LPT(1),'Q',x_LPT(1),fluid);
            rho_LPT(1) = CoolProp.PropsSI('D','T',T_LPT(1),'Q',x_LPT(1),fluid);
            u_LPT(1) = CoolProp.PropsSI('U','T',T_LPT(1),'Q',x_LPT(1),fluid);
            h_LPT(1) = CoolProp.PropsSI('H','T',T_LPT(1),'Q',x_LPT(1),fluid);
            s_LPT(1) = CoolProp.PropsSI('S','T',T_LPT(1),'Q',x_LPT(1),fluid);
        end
    end
    ex_LPT(1) = (u_LPT(1)-u_o) - T_o*(s_LPT(1)-s_o) + P_o*(1./rho_LPT(1)-1/rho_o);
    
%-------------------- Pump isoentropic efficiency ------------------------%
    if any(strcmp('eta_p',varargin))
        if ~isnumeric(varargin{find(strcmp('eta_p',varargin))+1})
            error('eta_p must be a number, not a %s.',class(varargin{find(strcmp('eta_p',varargin))+1}))
        end
        eta_p = varargin{find(strcmp('eta_p',varargin))+1};
    end

%--------------------- PRESSURE DROP IN THE HEATER -----------------------%
    if any(strcmp('DP_Hea',varargin))
        DeltaP_Hea = varargin{find(strcmp('DP_Hea',varargin))+1};
        if ~isnumeric(varargin{find(strcmp('DP_Hea',varargin))+1})
            error('DP_Hea must be a number, not a %s.',class(varargin{find(strcmp('DP_Hea',varargin))+1}))
        end
% ACRESCENTAR FORMAS DE CALCULAR A PERDA DE PRESSÃO EM FUNÇÃO DE OUTROS PARÂMETROS                        
    else
        DeltaP_Hea = 0;
    end    
%=========================== Discharging phase ===========================%
% States:                                                                 %
% 1 -> HPT outlet / evaporator inlet                                      %
% 2 -> evaporator outlet / turbine inlet                                  %
% 3 -> turbine outlet / condenser inlet                                   %
% 4 -> condenser outlet / LPT inlet                                       %
%-------------------------------------------------------------------------%

%------------------------ Parameters at t = 0 [s] ------------------------%
% Properties at the HPT output
    P_1(1) = P_HPT(1);
    if ph_HPT(1)==0 
        T_1(1) = T_HPT(1);
        h_1(1) = CoolProp.PropsSI('H','P',P_1(1),'T',T_1(1),fluid);
        s_1(1) = CoolProp.PropsSI('S','P',P_1(1),'T',T_1(1),fluid);
    elseif ph_HPT(1)==6
        T_1(1) = CoolProp.PropsSI('T','P',P_1(1),'Q',0,fluid);
        h_1(1) = CoolProp.PropsSI('H','P',P_1(1),'Q',0,fluid);
        s_1(1) = CoolProp.PropsSI('S','P',P_1(1),'Q',0,fluid);
    else
        error('The fluid at the HPT must be either liquid or a liquid-vapour mixture.')
    end
    
% Properties at the evaporator outlet / turbine inlet
    P_2(1) = P_1(1) - DeltaP_Ev;
    h_2(1) = CoolProp.PropsSI('H','P',P_2(1),'T',T_2(1),fluid);
    s_2(1) = CoolProp.PropsSI('S','P',P_2(1),'T',T_2(1),fluid);

% State at turbine outlet / condenser inlet
    P_3(1) = P_LPT(1);
    h_t(1) = CoolProp.PropsSI('H','P',P_3(1),'S',s_2(1),fluid);
    h_3(1) = turbine('h_in',h_2(1),'h_s',h_t(1),'eta_t',eta_t);
    T_3(1) = CoolProp.PropsSI('T','P',P_3(1),'H',h_3(1),fluid);
    s_3(1) = CoolProp.PropsSI('S','P',P_3(1),'H',h_3(1),fluid);
    
% State at condenser outlet / pump inlet
    P_4(1) = P_LPT(1); %!!!
    x_4(1) = 0;
    T_4(1) = CoolProp.PropsSI('T','P',P_4(1),'Q',x_4(1),fluid);
    h_4(1) = CoolProp.PropsSI('H','P',P_4(1),'Q',x_4(1),fluid);
    s_4(1) = CoolProp.PropsSI('S','P',P_4(1),'Q',x_4(1),fluid);    
    
% Energy transfers during the first time step    
    q_ev(1) = h_2(1) - h_1(1);
    w_t(1) = h_2(1) - h_3(1);

    switch Heat_source
        case 'HPs' % System with heat pump to provide heat at the evaporator/superheater
            % Hot source temperature (ORES system side of the heat exchanger)
            if T_2(1)>T_HPT(1) % If there is superheating
                h_g_ev(1) = CoolProp.PropsSI('H','P',P_2(1),'Q',1,fluid);
                q_SH(1) = h_2(1) - h_g_ev(1);
                q_lg(1) = q_ev(1) - q_SH(1);
                T_H1(1) = (q_lg(1)/q_ev(1))*T_HPT(1) + (q_SH(1)/q_ev(1))*(T_HPT(1)+T_2(1))/2;
            elseif T_2(1)==T_HPT(1) % If saturated state
                T_H1(1) = T_HPT(1);
            else
                error('T_2 < T_sat')
            end

            COP_carnot_HP1(1) = 1/(1-T_amb/T_H1(1));    
            if COP_HP1_0 == 0
                COP_HP1(1) = COP_carnot_HP1(1);
            elseif COP_HP1_0 < COP_carnot_HP1(1)
                COP_HP1(1) = COP_HP1_0;
            else
                error('COP_HP1 must be smaller than COP_carnot.')
            end
            w_c_HP1(1) = q_ev(1)/COP_HP1(1); % Work for reversible heat pump
            q_L_HP1(1) = q_ev(1) - w_c_HP1(1);      % Heat extracted from heat source
%             w_g(1) = w_t(1)-w_c_HP1(1);
            w_g(1) = w_t(1);
        case 'Free'
            w_g(1) = w_t(1);
    end

% Mass flow rate
    m_disch(1) = W_g/w_g(1);
    
    dE_t(1) = m_disch(1)*w_t(1)*dt;
    dQ_ev(1) = m_disch(1)*q_ev(1)*dt;

    Delta_m_disch_min = m_disch(1)*Delta_t;
    rho_HPT_end = CoolProp.PropsSI('D','P',P_HPT(1),'Q',x_HPT_f_disch,fluid);
    rho_LPT_end = CoolProp.PropsSI('D','P',P_LPT(1),'Q',x_LPT_f_disch,fluid);
    
    %!!!!!!!!!!!!!!!!!!!!!!! EM CONSTRUÇÃO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%---------------------- VOLUME OF THE STORAGE TANKS ----------------------%
    if any(strcmp('V_HPT',varargin)) % High pressure tank
        if ~isnumeric(varargin{find(strcmp('V_HPT',varargin))+1})
            error('V_HPT must be a number, not a %s.',class(varargin{find(strcmp('V_HPT',varargin))+1}))
        end
        V_HPT = varargin{find(strcmp('V_HPT',varargin))+1};
        % CHECAR CAPACIDADE DE MASSA
    elseif any(strcmp('K_V_HPT',varargin))    % INPUT K_V_HPT
        if ~isnumeric(varargin{find(strcmp('K_V_HPT',varargin))+1})
            error('K_V_HPT must be a number, not a %s.',class(varargin{find(strcmp('K_V_HPT',varargin))+1}))
        end
        K_V_H = varargin{find(strcmp('K_V_HPT',varargin))+1};
        V_HPT = K_V_H*( Delta_m_disch_min/(rho_HPT(1) - rho_HPT_end) );
    else
        V_HPT = K_V_H*( Delta_m_disch_min/(rho_HPT(1) - rho_HPT_end) );
    end 
        
    if any(strcmp('V_LPT',varargin)) % Low pressure tank
        if ~isnumeric(varargin{find(strcmp('V_LPT',varargin))+1})
            error('V_LPT must be a number, not a %s.',class(varargin{find(strcmp('V_LPT',varargin))+1}))
        end
        V_LPT = varargin{find(strcmp('V_LPT',varargin))+1};
    elseif any(strcmp('K_V_LPT',varargin))    % 
        if ~isnumeric(varargin{find(strcmp('K_V_LPT',varargin))+1})
            error('K_V_LPT must be a number, not a %s.',class(varargin{find(strcmp('K_V_LPT',varargin))+1}))
        end
        K_V_L = varargin{find(strcmp('K_V_LPT',varargin))+1};
%         V_LPT = K_V_L*V_HPT;
        V_LPT = K_V_L*( Delta_m_disch_min/(rho_LPT_end - rho_LPT(1)) );
%         CHECAR CAPACIDADE DE MASSA
    else
        V_LPT = K_V_L*V_HPT;
    end
    
    %!!!!!!!!!!!!!!!!!!!!!!! EM CONSTRUÇÃO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
%--------------------------------------------------------------------------    
    
% Total mass and internal energy at the HPT and LPT
    m_HPT(1) = V_HPT*rho_HPT(1);
    U_HPT(1) = m_HPT(1)*u_HPT(1);
    EX_HPT(1) = m_HPT(1)*ex_HPT(1);
    
    m_LPT(1) = V_LPT*rho_LPT(1);
    U_LPT(1) = m_LPT(1)*u_LPT(1);
    EX_LPT(1) = m_LPT(1)*ex_LPT(1);
    
    dE_net(1) = m_disch(1)*(w_t(1)-w_c_HP1(1))*dt;
    dQ_L_HP1(1) = m_disch(1)*q_L_HP1(1)*dt;
   
%-------------------- Configurações para plot dinâmico -------------------%
%     if plotBool
%         dyn = figure('color',[1 1 1]);
%         subplot(2,2,1)
%         hold on
%         subplot(2,2,2)
%         hold on
%         subplot(2,2,3)
%         hold on
%     end
%--------------------------- Dynamic analysis ----------------------------%
exc_count = 0;
exc_count2 = 0;
    for i=2:N_Steps+1
        t_disch(i) = t_disch(i-1) + dt;
    % HPT and LPT state update
        m_HPT(i) = m_HPT(i-1) - m_disch(i-1)*dt;
        rho_HPT(i) = m_HPT(i)/V_HPT;
        U_HPT(i) = U_HPT(i-1) - m_disch(i-1)*h_1(i-1)*dt;
        u_HPT(i) = U_HPT(i)/m_HPT(i);
        
        m_LPT(i) = m_LPT(i-1) + m_disch(i-1)*dt;
        rho_LPT(i) = m_LPT(i)/V_LPT;
        U_LPT(i) = U_LPT(i-1) + m_disch(i-1)*h_4(i-1)*dt;
        u_LPT(i) = U_LPT(i)/m_LPT(i);
        
        try
            if (rho_HPT(i)<rho_lim(1) || rho_HPT(i)>rho_lim(2)) && (rho_LPT(i)<rho_lim(1) || rho_LPT(i)>rho_lim(2))
                % If both tanks are outside unstable region
                P_HPT(i) = CoolProp.PropsSI('P','U',u_HPT(i),'D',rho_HPT(i),fluid);
                T_HPT(i) = CoolProp.PropsSI('T','U',u_HPT(i),'D',rho_HPT(i),fluid);
                x_HPT(i) = CoolProp.PropsSI('Q','U',u_HPT(i),'D',rho_HPT(i),fluid);
                h_HPT(i) = CoolProp.PropsSI('H','U',u_HPT(i),'D',rho_HPT(i),fluid);
                s_HPT(i) = CoolProp.PropsSI('S','U',u_HPT(i),'D',rho_HPT(i),fluid);
                ex_HPT(i) = (u_HPT(i)-u_o) - T_o*(s_HPT(i)-s_o) + P_o*(1./rho_HPT(i)-1/rho_o);
                EX_HPT(i) = m_HPT(i)*ex_HPT(i);

                P_LPT(i) = CoolProp.PropsSI('P','U',u_LPT(i),'D',rho_LPT(i),fluid);
                T_LPT(i) = CoolProp.PropsSI('T','U',u_LPT(i),'D',rho_LPT(i),fluid);
                x_LPT(i) = CoolProp.PropsSI('Q','U',u_LPT(i),'D',rho_LPT(i),fluid);
                h_LPT(i) = CoolProp.PropsSI('H','U',u_LPT(i),'D',rho_LPT(i),fluid);
                s_LPT(i) = CoolProp.PropsSI('S','U',u_LPT(i),'D',rho_LPT(i),fluid);
                ex_LPT(i) = (u_LPT(i)-u_o) - T_o*(s_LPT(i)-s_o) + P_o*(1./rho_LPT(i)-1/rho_o);
                EX_LPT(i) = m_LPT(i)*ex_LPT(i);

            elseif (rho_HPT(i)>rho_lim(1) && rho_HPT(i)<rho_lim(2)) && (rho_LPT(i)<rho_lim(1) || rho_LPT(i)>rho_lim(2))
                % If only the HPT is inside unstable region
                if (rho_HPT(i)>rho_lim(1) && rho_HPT(i)<rho_lim(2)) && (rho_HPT(i-1)<rho_lim(1) || rho_HPT(i-1)>rho_lim(2))
                    % If HPT was outside the unstable region in the last iteration
    %                 disp('update_HPT')
                    % Set the superior limit state
                    u_lim_s_HPT = u_HPT(i-1);
                    rho_lim_s_HPT = rho_HPT(i-1);
                    P_lim_s_HPT = CoolProp.PropsSI('P','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
                    T_lim_s_HPT = CoolProp.PropsSI('T','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
                    x_lim_s_HPT = CoolProp.PropsSI('Q','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
                    h_lim_s_HPT = CoolProp.PropsSI('H','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
                    s_lim_s_HPT = CoolProp.PropsSI('S','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);

                    % Calculate the required mass to exit the unstable region
                    dm = (rho_HPT(i-1) - rho_lim(1))*V_HPT;
                    m_lim_i_HPT = m_HPT(i-1) - dm;

                    % Set the inferior limit state
                    rho_lim_i_HPT = m_lim_i_HPT/V_HPT;
                    U_lim_i_HPT = U_HPT(i-1) - dm*h_1(i-1);
                    u_lim_i_HPT = U_lim_i_HPT/m_lim_i_HPT;
                    P_lim_i_HPT = CoolProp.PropsSI('P','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                    T_lim_i_HPT = CoolProp.PropsSI('T','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                    x_lim_i_HPT = CoolProp.PropsSI('Q','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                    h_lim_i_HPT = CoolProp.PropsSI('H','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                    s_lim_i_HPT = CoolProp.PropsSI('S','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                end

                % Linearizarion of the unstable region for the HPT
                P_HPT(i) = P_lim_i_HPT + (rho_HPT(i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(P_lim_s_HPT - P_lim_i_HPT);
                T_HPT(i) = T_lim_i_HPT + (rho_HPT(i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(T_lim_s_HPT - T_lim_i_HPT);
                x_HPT(i) = x_lim_i_HPT + (rho_HPT(i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(x_lim_s_HPT - x_lim_i_HPT);
                h_HPT(i) = h_lim_i_HPT + (rho_HPT(i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(h_lim_s_HPT - h_lim_i_HPT);
                s_HPT(i) = s_lim_i_HPT + (rho_HPT(i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(s_lim_s_HPT - s_lim_i_HPT);
                ex_HPT(i) = (u_HPT(i)-u_o) - T_o*(s_HPT(i)-s_o) + P_o*(1./rho_HPT(i)-1/rho_o);
                EX_HPT(i) = m_HPT(i)*ex_HPT(i);

                P_LPT(i) = CoolProp.PropsSI('P','U',u_LPT(i),'D',rho_LPT(i),fluid);
                T_LPT(i) = CoolProp.PropsSI('T','U',u_LPT(i),'D',rho_LPT(i),fluid);
                x_LPT(i) = CoolProp.PropsSI('Q','U',u_LPT(i),'D',rho_LPT(i),fluid);
                h_LPT(i) = CoolProp.PropsSI('H','U',u_LPT(i),'D',rho_LPT(i),fluid);
                s_LPT(i) = CoolProp.PropsSI('S','U',u_LPT(i),'D',rho_LPT(i),fluid);
                ex_LPT(i) = (u_LPT(i)-u_o) - T_o*(s_LPT(i)-s_o) + P_o*(1./rho_LPT(i)-1/rho_o);
                EX_LPT(i) = m_LPT(i)*ex_LPT(i);
            elseif (rho_HPT(i)<rho_lim(1) || rho_HPT(i)>rho_lim(2)) && (rho_LPT(i)>rho_lim(1) && rho_LPT(i)<rho_lim(2))
                % If only the LPT is inside unstable region
                P_HPT(i) = CoolProp.PropsSI('P','U',u_HPT(i),'D',rho_HPT(i),fluid);
                T_HPT(i) = CoolProp.PropsSI('T','U',u_HPT(i),'D',rho_HPT(i),fluid);
                x_HPT(i) = CoolProp.PropsSI('Q','U',u_HPT(i),'D',rho_HPT(i),fluid);
                h_HPT(i) = CoolProp.PropsSI('H','U',u_HPT(i),'D',rho_HPT(i),fluid);
                s_HPT(i) = CoolProp.PropsSI('S','U',u_HPT(i),'D',rho_HPT(i),fluid);
                ex_HPT(i) = (u_HPT(i)-u_o) - T_o*(s_HPT(i)-s_o) + P_o*(1./rho_HPT(i)-1/rho_o);
                EX_HPT(i) = m_HPT(i)*ex_HPT(i);

                if (rho_LPT(i)>rho_lim(1) && rho_LPT(i)<rho_lim(2)) && (rho_LPT(i-1)<rho_lim(1) || rho_LPT(i-1)>rho_lim(2))
                    % If LPT was outside the unstable region in the last iteration
    %                 disp('update_LPT')
                    % Set the inferior limit state
                    u_lim_i_LPT = u_LPT(i-1);
                    rho_lim_i_LPT = rho_LPT(i-1);
                    P_lim_i_LPT = CoolProp.PropsSI('P','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
                    T_lim_i_LPT = CoolProp.PropsSI('T','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
                    x_lim_i_LPT = CoolProp.PropsSI('Q','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
                    h_lim_i_LPT = CoolProp.PropsSI('H','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
                    s_lim_i_LPT = CoolProp.PropsSI('S','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);

                    % Estimate required mass to exit unstable region
                    dm = (rho_lim(2) - rho_LPT(i-1))*V_LPT;
                    m_lim_s_LPT = m_LPT(i-1) + dm;

                    % Set the superior limit state
                    rho_lim_s_LPT = m_lim_s_LPT/V_LPT;
                    U_lim_s_LPT = U_LPT(i-1) + dm*h_4(i-1);
                    u_lim_s_LPT = U_lim_s_LPT/m_lim_s_LPT;
                    P_lim_s_LPT = CoolProp.PropsSI('P','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                    T_lim_s_LPT = CoolProp.PropsSI('T','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                    x_lim_s_LPT = CoolProp.PropsSI('Q','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                    h_lim_s_LPT = CoolProp.PropsSI('H','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                    s_lim_s_LPT = CoolProp.PropsSI('S','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                end

                % Linearizarion of the unstable region for the LPT
                P_LPT(i) = P_lim_i_LPT + (rho_LPT(i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(P_lim_s_LPT - P_lim_i_LPT);
                T_LPT(i) = T_lim_i_LPT + (rho_LPT(i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(T_lim_s_LPT - T_lim_i_LPT);
                x_LPT(i) = x_lim_i_LPT + (rho_LPT(i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(x_lim_s_LPT - x_lim_i_LPT);
                h_LPT(i) = h_lim_i_LPT + (rho_LPT(i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(h_lim_s_LPT - h_lim_i_LPT);
                s_LPT(i) = s_lim_i_LPT + (rho_LPT(i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(s_lim_s_LPT - s_lim_i_LPT);
                ex_LPT(i) = (u_LPT(i)-u_o) - T_o*(s_LPT(i)-s_o) + P_o*(1./rho_LPT(i)-1/rho_o);
                EX_LPT(i) = m_LPT(i)*ex_LPT(i);

            elseif (rho_HPT(i)>rho_lim(1) && rho_HPT(i)<rho_lim(2)) && (rho_LPT(i)>rho_lim(1) && rho_LPT(i)<rho_lim(2))
                % If both tanks are in the unstable region
                if (rho_HPT(i)>rho_lim(1) && rho_HPT(i)<rho_lim(2)) && (rho_HPT(i-1)<rho_lim(1) || rho_HPT(i-1)>rho_lim(2))
                    % If the HPT was outside the unstable region in the last iteration 
    %                 disp('update_HPT 2')

                    % Set the superior limit state
                    u_lim_s_HPT = u_HPT(i-1);
                    rho_lim_s_HPT = rho_HPT(i-1);
                    P_lim_s_HPT = CoolProp.PropsSI('P','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
                    T_lim_s_HPT = CoolProp.PropsSI('T','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
                    x_lim_s_HPT = CoolProp.PropsSI('Q','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
                    h_lim_s_HPT = CoolProp.PropsSI('H','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
                    s_lim_s_HPT = CoolProp.PropsSI('S','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);

                    % Estimate required mass to exit unstable region
                    dm = (rho_HPT(i-1) - rho_lim(1))*V_HPT;
                    m_lim_i_HPT = m_HPT(i-1) - dm;

                    % Set the inferior limit state
                    rho_lim_i_HPT = m_lim_i_HPT/V_HPT;
                    U_lim_i_HPT = U_HPT(i-1) - dm*h_1(i-1);
                    u_lim_i_HPT = U_lim_i_HPT/m_lim_i_HPT;
                    P_lim_i_HPT = CoolProp.PropsSI('P','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                    T_lim_i_HPT = CoolProp.PropsSI('T','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                    x_lim_i_HPT = CoolProp.PropsSI('Q','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                    h_lim_i_HPT = CoolProp.PropsSI('H','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                    s_lim_i_HPT = CoolProp.PropsSI('S','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                end

                % Linearizarion of the unstable region for the HPT
                P_HPT(i) = P_lim_i_HPT + (rho_HPT(i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(P_lim_s_HPT - P_lim_i_HPT);
                T_HPT(i) = T_lim_i_HPT + (rho_HPT(i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(T_lim_s_HPT - T_lim_i_HPT);
                x_HPT(i) = x_lim_i_HPT + (rho_HPT(i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(x_lim_s_HPT - x_lim_i_HPT);
                h_HPT(i) = h_lim_i_HPT + (rho_HPT(i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(h_lim_s_HPT - h_lim_i_HPT);
                s_HPT(i) = s_lim_i_HPT + (rho_HPT(i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(s_lim_s_HPT - s_lim_i_HPT);
                ex_HPT(i) = (u_HPT(i)-u_o) - T_o*(s_HPT(i)-s_o) + P_o*(1./rho_HPT(i)-1/rho_o);
                EX_HPT(i) = m_HPT(i)*ex_HPT(i);

                if (rho_LPT(i)>rho_lim(1) && rho_LPT(i)<rho_lim(2)) && (rho_LPT(i-1)<rho_lim(1) || rho_LPT(i-1)>rho_lim(2))
                    % If the LPT was outside the unstable region in the last iteration 
    %                 disp('update_LPT 2')

                    % Set the inferior limit state
                    u_lim_i_LPT = u_LPT(i-1);
                    rho_lim_i_LPT = rho_LPT(i-1);
                    P_lim_i_LPT = CoolProp.PropsSI('P','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
                    T_lim_i_LPT = CoolProp.PropsSI('T','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
                    x_lim_i_LPT = CoolProp.PropsSI('Q','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
                    h_lim_i_LPT = CoolProp.PropsSI('H','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
                    s_lim_i_LPT = CoolProp.PropsSI('S','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);

                    % Estimate required mass to exit unstable region
                    dm = (rho_lim(2) - rho_LPT(i-1))*V_LPT;
                    m_lim_s_LPT = m_LPT(i-1) + dm;

                    % Set the superior limit state
                    rho_lim_s_LPT = m_lim_s_LPT/V_LPT;
                    U_lim_s_LPT = U_LPT(i-1) + dm*h_4(i-1);
                    u_lim_s_LPT = U_lim_s_LPT/m_lim_s_LPT;
                    P_lim_s_LPT = CoolProp.PropsSI('P','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                    T_lim_s_LPT = CoolProp.PropsSI('T','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                    x_lim_s_LPT = CoolProp.PropsSI('Q','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                    h_lim_s_LPT = CoolProp.PropsSI('H','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                    s_lim_s_LPT = CoolProp.PropsSI('S','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                end

                % Linearizarion of the unstable region for the LPT
                P_LPT(i) = P_lim_i_LPT + (rho_LPT(i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(P_lim_s_LPT - P_lim_i_LPT);
                T_LPT(i) = T_lim_i_LPT + (rho_LPT(i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(T_lim_s_LPT - T_lim_i_LPT);
                x_LPT(i) = x_lim_i_LPT + (rho_LPT(i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(x_lim_s_LPT - x_lim_i_LPT);
                h_LPT(i) = h_lim_i_LPT + (rho_LPT(i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(h_lim_s_LPT - h_lim_i_LPT);
                s_LPT(i) = s_lim_i_LPT + (rho_LPT(i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(s_lim_s_LPT - s_lim_i_LPT);
                ex_LPT(i) = (u_LPT(i)-u_o) - T_o*(s_LPT(i)-s_o) + P_o*(1./rho_LPT(i)-1/rho_o);
                EX_LPT(i) = m_LPT(i)*ex_LPT(i);
            end
            
                        
            
            exc_count = 0;
        catch ME
            switch ME.identifier
                case 'SWIG:RuntimeError'
                    if i>2
                    
                        P_HPT(i) = P_HPT(i-1) + (rho_HPT(i)-rho_HPT(i-1))/(rho_HPT(i-1)-rho_HPT(i-2))*(P_HPT(i-1) - P_HPT(i-2));
                        T_HPT(i) = T_HPT(i-1) + (rho_HPT(i)-rho_HPT(i-1))/(rho_HPT(i-1)-rho_HPT(i-2))*(T_HPT(i-1) - T_HPT(i-2));
                        x_HPT(i) = x_HPT(i-1) + (rho_HPT(i)-rho_HPT(i-1))/(rho_HPT(i-1)-rho_HPT(i-2))*(x_HPT(i-1) - x_HPT(i-2));
                        h_HPT(i) = h_HPT(i-1) + (rho_HPT(i)-rho_HPT(i-1))/(rho_HPT(i-1)-rho_HPT(i-2))*(h_HPT(i-1) - h_HPT(i-2));
                        s_HPT(i) = s_HPT(i-1) + (rho_HPT(i)-rho_HPT(i-1))/(rho_HPT(i-1)-rho_HPT(i-2))*(s_HPT(i-1) - s_HPT(i-2));

                        ex_HPT(i) = (u_HPT(i)-u_o) - T_o*(s_HPT(i)-s_o) + P_o*(1./rho_HPT(i)-1/rho_o);
                        EX_HPT(i) = m_HPT(i)*ex_HPT(i);

                        P_LPT(i) = P_LPT(i-1) + (rho_LPT(i)-rho_LPT(i-1))/(rho_LPT(i-1)-rho_LPT(i-2))*(P_LPT(i-1) - P_LPT(i-2));
                        T_LPT(i) = T_LPT(i-1) + (rho_LPT(i)-rho_LPT(i-1))/(rho_LPT(i-1)-rho_LPT(i-2))*(T_LPT(i-1) - T_LPT(i-2));
                        x_LPT(i) = x_LPT(i-1) + (rho_LPT(i)-rho_LPT(i-1))/(rho_LPT(i-1)-rho_LPT(i-2))*(x_LPT(i-1) - x_LPT(i-2));
                        h_LPT(i) = h_LPT(i-1) + (rho_LPT(i)-rho_LPT(i-1))/(rho_LPT(i-1)-rho_LPT(i-2))*(h_LPT(i-1) - h_LPT(i-2));
                        s_LPT(i) = s_LPT(i-1) + (rho_LPT(i)-rho_LPT(i-1))/(rho_LPT(i-1)-rho_LPT(i-2))*(s_LPT(i-1) - s_LPT(i-2));

                        ex_LPT(i) = (u_LPT(i)-u_o) - T_o*(s_LPT(i)-s_o) + P_o*(1./rho_LPT(i)-1/rho_o);
                        EX_LPT(i) = m_LPT(i)*ex_LPT(i);
                    else
                        P_HPT(i) = P_HPT(i-1);
                        T_HPT(i) = T_HPT(i-1);
                        x_HPT(i) = x_HPT(i-1);
                        h_HPT(i) = h_HPT(i-1);
                        s_HPT(i) = s_HPT(i-1);

                        ex_HPT(i) = (u_HPT(i)-u_o) - T_o*(s_HPT(i)-s_o) + P_o*(1./rho_HPT(i)-1/rho_o);
                        EX_HPT(i) = m_HPT(i)*ex_HPT(i);

                        P_LPT(i) = P_LPT(i-1);
                        T_LPT(i) = T_LPT(i-1);
                        x_LPT(i) = x_LPT(i-1);
                        h_LPT(i) = h_LPT(i-1);
                        s_LPT(i) = s_LPT(i-1);

                        ex_LPT(i) = (u_LPT(i)-u_o) - T_o*(s_LPT(i)-s_o) + P_o*(1./rho_LPT(i)-1/rho_o);
                        EX_LPT(i) = m_LPT(i)*ex_LPT(i);
                    end
                    exc_count = exc_count + 1;
                    if exc_count >= 5
                       error('5 consecutive errors.') 
                    end
                otherwise
                    rethrow(ME)
            end
        end         
        
    % Cycle state update
    % State at HPT outlet / Evaporator inlet
        P_1(i) = P_HPT(i);
%         if ph_HPT(i)==0
%             T_1(i) = T_HPT(i);
%             h_1(i) = CoolProp.PropsSI('H','P',P_1(i),'T',T_1(i),fluid);
%             s_1(i) = CoolProp.PropsSI('S','P',P_1(i),'T',T_1(i),fluid);
%         elseif ph_HPT(i)==6
            T_1(i) = CoolProp.PropsSI('T','P',P_1(i),'Q',0,fluid);
            h_1(i) = CoolProp.PropsSI('H','P',P_1(i),'Q',0,fluid);
            s_1(i) = CoolProp.PropsSI('S','P',P_1(i),'Q',0,fluid);
%         else
%             error('The fluid at the HPT must be either liquid or a liquid-vapour mixture.')
%         end
    % State at Evaporator outlet / Turbine inlet
        P_2(i) = P_1(i) - DeltaP_Ev;
        T_sat_Ev(i) = CoolProp.PropsSI('T','P',P_2(i),'Q',0,fluid);
        T_2(i) = T_sat_Ev(i) + DT_SH;
        h_2(i) = CoolProp.PropsSI('H','P',P_2(i),'T',T_2(i),fluid);
        s_2(i) = CoolProp.PropsSI('S','P',P_2(i),'T',T_2(i),fluid);
    
    % State at turbine outlet / condenser inlet
        P_3(i) = P_LPT(i);
        h_t(i) = CoolProp.PropsSI('H','P',P_3(i),'S',s_2(i),fluid);
        h_3(i) = turbine('h_in',h_2(i),'h_s',h_t(i),'eta_t',eta_t);
        T_3(i) = CoolProp.PropsSI('T','P',P_3(i),'H',h_3(i),fluid);
        s_3(i) = CoolProp.PropsSI('S','P',P_3(i),'H',h_3(i),fluid);
        
    % State at condenser outlet / LPT inlet
        P_4(i) = P_LPT(i); %!!!
        x_4(i) = 0;
        T_4(i) = CoolProp.PropsSI('T','P',P_4(i),'Q',x_4(i),fluid);
        h_4(i) = CoolProp.PropsSI('H','P',P_4(i),'Q',x_4(i),fluid);
        s_4(i) = CoolProp.PropsSI('S','P',P_4(i),'Q',x_4(i),fluid);    
    
    % Energy transfers during time step i
        q_ev(i) = h_2(i) - h_1(i);
        %q_he(1) = h_1(1) - h_5(1);
        %q_in(1) = q_he(1) + q_ev(1);
        w_t(i) = h_2(i) - h_3(i);
        
%         if i == 474
%             i
%         end
        switch Heat_source
            case 'HPs' % System with heat pump to provide heat at the evaporator/superheater
                % Hot source temperature (ORES system side of the heat exchanger)
                if T_2(i)>T_HPT(i) % If there is superheating
                    h_g_ev(i) = CoolProp.PropsSI('H','P',P_2(i),'Q',1,fluid);
                    q_SH(i) = h_2(i) - h_g_ev(i);
                    q_lg(i) = q_ev(i) - q_SH(i);
                    T_H1(i) = (q_lg(i)/q_ev(i))*T_HPT(i) + (q_SH(i)/q_ev(i))*(T_HPT(i)+T_2(i))/2;
                    exc_count2 = 0;
                elseif T_2(i)==T_HPT(i) % If saturated state
                    T_H1(i) = T_HPT(i);
                    exc_count2 = 0;
                else
                    exc_count2 = exc_count2 + 1;
                    warning('T_2 < T_sat')
                    T_H1(i) = T_H1(i-1) + (T_2(i)-T_2(i-1))/(T_2(i-1)-T_2(i-2))*(T_H1(i-1) - T_H1(i-2));
%                     [T_2(i) T_sat_Ev(i)]
                    if exc_count2 > 5
                        error('T_2 < T_sat')
                    end
                end

                COP_carnot_HP1(i) = 1/(1-T_amb/T_H1(i));
                if COP_HP1_0 == 0
                    COP_HP1(i) = COP_carnot_HP1(i);
                elseif COP_HP1_0 < COP_carnot_HP1(i)
                    COP_HP1(i) = COP_HP1_0;
                else
                    error('COP_HP1 must be smaller than COP_carnot.')
                end
                w_c_HP1(i) = q_ev(i)/COP_HP1(i); % Work for reversible heat pump
                q_L_HP1(i) = q_ev(i) - w_c_HP1(i);      % Heat extracted from heat source
%                 w_g(i) = w_t(i)-w_c_HP1(i);
                w_g(i) = w_t(i);
            case 'Free'
                w_g(i) = w_t(i);
        end
              
        m_disch(i) = W_g/w_g(i);
        dE_t(i) = m_disch(i)*w_t(i)*dt;
        dQ_ev(i) = m_disch(i)*q_ev(i)*dt;

        dE_net(i) = m_disch(i)*(w_t(i)-w_c_HP1(i))*dt;
        dQ_L_HP1(i) = m_disch(i)*q_L_HP1(i)*dt;

        if x_LPT(i)<x_LPT_f_disch
            disp('x_LPT < 0.02 during discharge.')
            break;
        elseif x_HPT(i)>=x_HPT_f_disch
            disp('x_HPT >= 0.98 during discharge.')
            break;
        elseif P_HPT(i) < P_lim_inf
            disp('P_HPT < P_lim.')
            break;
        elseif w_t(i)<w_c_HP1(i)
            disp('Generated power < Consumed power.')
            break;
        end
    
    end
n_steps_disch = i;
t_generation = t_disch(i);
if t_generation < Delta_t
    t_disch = t_disch(1:n_steps_disch);
    P_LPT = P_LPT(1:n_steps_disch);
    x_LPT = x_LPT(1:n_steps_disch);
    T_LPT = T_LPT(1:n_steps_disch);
    rho_LPT = rho_LPT(1:n_steps_disch);
    u_LPT = u_LPT(1:n_steps_disch);
    h_LPT = h_LPT(1:n_steps_disch);
    s_LPT = s_LPT(1:n_steps_disch);
    ex_LPT = ex_LPT(1:n_steps_disch);
    m_LPT = m_LPT(1:n_steps_disch);
    U_LPT = U_LPT(1:n_steps_disch);
    EX_LPT = EX_LPT(1:n_steps_disch);
    
    P_HPT = P_HPT(1:n_steps_disch);
    x_HPT = x_HPT(1:n_steps_disch);
    T_HPT = T_HPT(1:n_steps_disch);
    rho_HPT = rho_HPT(1:n_steps_disch);
    u_HPT = u_HPT(1:n_steps_disch);
    h_HPT = h_HPT(1:n_steps_disch);
    s_HPT = s_HPT(1:n_steps_disch);
    ex_HPT = ex_HPT(1:n_steps_disch);
    m_HPT = m_HPT(1:n_steps_disch);
    U_HPT = U_HPT(1:n_steps_disch);
    EX_HPT = EX_HPT(1:n_steps_disch);
    
    P_1 = P_1(1:n_steps_disch);
    T_1 = T_1(1:n_steps_disch);
    s_1 = s_1(1:n_steps_disch);
    h_1 = h_1(1:n_steps_disch);
    
    P_2 = P_2(1:n_steps_disch);
    T_2 = T_2(1:n_steps_disch);
    s_2 = s_2(1:n_steps_disch);
    h_2 = h_2(1:n_steps_disch);
    
    P_3 = P_3(1:n_steps_disch);
    T_3 = T_3(1:n_steps_disch);
    s_3 = s_3(1:n_steps_disch);
    h_t = h_t(1:n_steps_disch);
    h_3 = h_3(1:n_steps_disch);
    
    P_4 = P_4(1:n_steps_disch);
    T_4 = T_4(1:n_steps_disch);
    s_4 = s_4(1:n_steps_disch);
    h_4 = h_4(1:n_steps_disch);
    x_4 = x_4(1:n_steps_disch);
    
    q_ev = q_ev(1:n_steps_disch);
    w_t = w_t(1:n_steps_disch);
    w_g = w_g(1:n_steps_disch);
            
%     dE_g = dE_g(1:n_steps_disch);
    dE_t = dE_t(1:n_steps_disch);
    dQ_ev = dQ_ev(1:n_steps_disch);
    
    dE_net = dE_net(1:n_steps_disch);

    switch Heat_source
        case 'HPs'
            h_g_ev = h_g_ev(1:n_steps_disch);
            q_SH = q_SH(1:n_steps_disch);
            q_lg = q_lg(1:n_steps_disch);
            T_H1 = T_H1(1:n_steps_disch);
            COP_carnot_HP1 = COP_carnot_HP1(1:n_steps_disch);
            COP_HP1 = COP_HP1(1:n_steps_disch);
            w_c_HP1 = w_c_HP1(1:n_steps_disch);
            q_L_HP1 = q_L_HP1(1:n_steps_disch);
            dQ_L_HP1 = dQ_L_HP1(1:n_steps_disch);
        otherwise
            h_g_ev = h_g_ev(1:n_steps_disch);
            q_SH = q_SH(1:n_steps_disch);
            q_lg = q_lg(1:n_steps_disch);
            T_H1 = T_H1(1:n_steps_disch);
            COP_carnot_HP1 = COP_carnot_HP1(1:n_steps_disch);
            COP_HP1 = COP_HP1(1:n_steps_disch);
            w_c_HP1 = w_c_HP1(1:n_steps_disch);
            q_L_HP1 = q_L_HP1(1:n_steps_disch);
            dQ_L_HP1 = dQ_L_HP1(1:n_steps_disch);
%             dE_HP1 = dE_HP1(1:n_steps_disch);
    end
    m_disch = m_disch(1:n_steps_disch);
    
    warning('t_generation smaller than desired Delta_t.')
    fprintf('t_generation = %d. \n',t_generation)
end
% P_HPT(end);
% P_LPT(end);

%---------------------------- Results display ----------------------------%
if plotBool
    % Properties over time
    figure('color',[1 1 1])
    %suptitle('V_HPT = ',V_HPT,'V_LPT = ',V_LPT)
    subplot(3,2,1)
    plot(t_disch,m_disch,'k')
    grid on;
    ylabel('Mass flow rate [kg/s]')
    grid on
    xtickformat('%,0.0f')
    ytickformat('%0.1f')
    %hl = title('$\dot{m} [kg-s^{-1}]$');
    %set(hl, 'Interpreter', 'latex');
    subplot(3,2,2)
    plot(t_disch,m_HPT(1:n_steps_disch),'k'); hold on;
    plot(t_disch,m_LPT(1:n_steps_disch),'color',[0.7 0.7 0.7]); 
    grid on;
    ylabel('Total mass [kg]')
    xtickformat('%,0.0f')
    ytickformat('%0.1f')
    legend('m_{HPT}','m_{LPT}','location','west')
    
    subplot(3,2,3)
    plot(t_disch,P_HPT(1:n_steps_disch)./1000,'k'); hold on;
    plot(t_disch,P_LPT(1:n_steps_disch)./1000,'color',[0.7 0.7 0.7]); 
    grid on;    
    ylabel('Pressure [kPa]')
    xtickformat('%,0.0f')
    ytickformat('%,0.0f')
    legend('P_{HPT}','P_{LPT}','location','west')
    
    subplot(3,2,4)
    plot(t_disch,T_HPT(1:n_steps_disch),'k'); hold on;
    plot(t_disch,T_LPT(1:n_steps_disch),'color',[0.7 0.7 0.7]); 
    grid on;
    ylabel('Temperature [K]')
    xtickformat('%,0.0f')
    legend('T_{HPT}','T_{LPT}','location','west')
    
    subplot(3,2,5)
    plot(t_disch,x_HPT(1:n_steps_disch),'k'); hold on;
    plot(t_disch,x_LPT(1:n_steps_disch),'color',[0.7 0.7 0.7]); 
    grid on;
    ylim([0 1])
    xtickformat('%,0.0f')
    ytickformat('%0.1f')
    legend('x_{HPT}','x_{LPT}','location','northeast')
    ylabel('Quality [-]')
    
    subplot(3,2,6)
    plot(t_disch,rho_HPT(1:n_steps_disch),'k'); hold on;
    plot(t_disch,rho_LPT(1:n_steps_disch),'color',[0.7 0.7 0.7]);
    grid on;
    ylabel('Specific mass [kg/m^3]')
    xlabel('Time [s]')
    xtickformat('%,0.0f')
    ytickformat('%0.0f')
    legend('\rho_{HPT}','\rho_{LPT}','location','west')
    
    % Plot the T-s diagramm for the fluid with the defined isobaric lines
%     TSDiag=plotTS(fluid,'Iso_P',[P_HPT(1),P_LPT(1),P_HPT(end),P_LPT(end)]);
    
    plotTS(fluid,'Iso_P',[P_HPT(1),P_LPT(1),P_HPT(n_steps_disch),P_LPT(n_steps_disch)]);
    
    plot(s_HPT(1:n_steps_disch)./1000,T_HPT(1:n_steps_disch),'k','LineWidth',2)
    handle_disch(1) =plot(s_LPT(1:n_steps_disch)./1000,T_LPT(1:n_steps_disch),'k','LineWidth',2);
    
    plotTS(fluid,'Process','Isobaric',P_1(1),s_1(1),s_2(1),'plotColor','b');%'k');
    plotTS(fluid,'Process','Expansion',P_2(1),h_2(1),P_3(1),eta_t,'plotColor','b');%,'k');
    handle_disch(2) = plotTS(fluid,'Process','Isobaric',P_3(1),s_3(1),s_4(1),'plotColor','b');%,'k');
    
    plotTS(fluid,'Process','Isobaric',P_1(n_steps_disch),s_1(n_steps_disch),s_2(n_steps_disch),'plotColor','r');%,[0.7 0.7 0.7]);
    plotTS(fluid,'Process','Expansion',P_1(n_steps_disch),h_2(n_steps_disch),P_3(n_steps_disch),eta_t,'plotColor','r');%,[0.7 0.7 0.7]);
    handle_disch(3) = plotTS(fluid,'Process','Isobaric',P_3(n_steps_disch),s_3(n_steps_disch),s_4(n_steps_disch),'plotColor','r');%,[0.7 0.7 0.7]);
    grid on;
    % T-s diagramm during discharging phase
    
    ylim([T_LPT(1)-10 T_2(1)+20]);
    xlim([floor(s_LPT(n_steps_disch)./1000) ceil(s_3(n_steps_disch)./100)./10])
    
    text([s_1(1)-50;s_2(1)+20;s_3(1)+20;s_4(1)-10]./1000,[T_1(1);T_2(1);T_3(1);T_4(1)+10],{'1','2','3','4'},'FontSize',12)
%     text([s_1(n_steps_disch);s_2(n_steps_disch);s_3(n_steps_disch);s_4(n_steps_disch)]./1000,...
%         [T_1(n_steps_disch);T_2(n_steps_disch);T_3(n_steps_disch);T_4(n_steps_disch)],{'1','2','3','4'},...
%         'Color',[0.5 0.5 0.5]);
    plot([s_1(1);s_2(1);s_3(1);s_4(1)]./1000,[T_1(1);T_2(1);T_3(1);T_4(1)],'b.','MarkerSize',20)
    plot([s_1(n_steps_disch);s_2(n_steps_disch);s_3(n_steps_disch);s_4(n_steps_disch)]./1000,...
         [T_1(n_steps_disch);T_2(n_steps_disch);T_3(n_steps_disch);T_4(n_steps_disch)],...
         'k.','MarkerSize',20,'Color','r')
    legend(handle_disch,{'State at the tanks','Process diagramm at t=0 s',...
        'Process diagramm at t=t_{disch}'},'location','northwest')
    xtickformat('%0.1f')
    
    fig = gcf;
    fig.PaperOrientation = 'landscape';
    fig.PaperSize = [6 4.5];
    
    pause(0.1)
end

%============================= Charging phase ============================%
% States:                                                                 %
% 5 -> LPT outlet / pump inlet                                            %
% 6 -> pump outlet / heater inlet                                         %
% 7 -> heater outlet / HPT inlet                                          %
%-------------------------------------------------------------------------%

%--------------------- Charging simulation parameters --------------------%
Delta_m_disch = m_HPT(1) - m_HPT(n_steps_disch);

%------------------------ Parameters at t = 0 [s] ------------------------%
% HPT and LPT state update
m_HPT(n_steps_disch+1) = m_HPT(n_steps_disch);
rho_HPT(n_steps_disch+1) = m_HPT(n_steps_disch+1)/V_HPT;
U_HPT(n_steps_disch+1) = U_HPT(n_steps_disch);
u_HPT(n_steps_disch+1) = u_HPT(n_steps_disch);
P_HPT(n_steps_disch+1) = P_HPT(n_steps_disch);
T_HPT(n_steps_disch+1) = T_HPT(n_steps_disch);
x_HPT(n_steps_disch+1) = x_HPT(n_steps_disch);
h_HPT(n_steps_disch+1) = h_HPT(n_steps_disch);
s_HPT(n_steps_disch+1) = s_HPT(n_steps_disch);
ex_HPT(n_steps_disch+1) = ex_HPT(n_steps_disch);
EX_HPT(n_steps_disch+1) = EX_HPT(n_steps_disch);

m_LPT(n_steps_disch+1) = m_LPT(n_steps_disch);
rho_LPT(n_steps_disch+1) = rho_LPT(n_steps_disch);
U_LPT(n_steps_disch+1) = U_LPT(n_steps_disch);
u_LPT(n_steps_disch+1) = u_LPT(n_steps_disch);
P_LPT(n_steps_disch+1) = P_LPT(n_steps_disch);
T_LPT(n_steps_disch+1) = T_LPT(n_steps_disch);
x_LPT(n_steps_disch+1) = x_LPT(n_steps_disch);
h_LPT(n_steps_disch+1) = h_LPT(n_steps_disch);
s_LPT(n_steps_disch+1) = s_LPT(n_steps_disch);
ex_LPT(n_steps_disch+1) = ex_LPT(n_steps_disch);
EX_LPT(n_steps_disch+1) = EX_LPT(n_steps_disch);

% m_HPT_0_c = m_HPT(n_steps_disch+1);

t_charg(1) = 0;

% State at LPT outlet / Pump inlet
P_5(1) = P_LPT(n_steps_disch); % INCLUIR VERIFICAÇÃO DE ESTADO LIQUIDO SATURADO OU SUBRESFRIADO
x_5(1) = 0;
T_5(1) = CoolProp.PropsSI('T','P',P_5(1),'Q',x_5(1),fluid);
h_5(1) = CoolProp.PropsSI('H','P',P_5(1),'Q',x_5(1),fluid);
s_5(1) = CoolProp.PropsSI('S','P',P_5(1),'Q',x_5(1),fluid);

% State at Pump outlet / Heater inlet
% P_6(1) = P_HPT(1); % VERIFICAR
P_6(1) = P_HPT(n_steps_disch);
h_p(1) = CoolProp.PropsSI('H','P',P_6(1),'S',s_5(1),fluid);
h_6(1) = pump('h_in',h_5(1),'h_s',h_p(1),'eta_p',eta_p);
T_6(1) = CoolProp.PropsSI('T','P',P_6(1),'H',h_6(1),fluid);
s_6(1) = CoolProp.PropsSI('S','P',P_6(1),'H',h_6(1),fluid);    

% State at Heater outlet / HPT inlet
P_7(1) = P_6(1) - DeltaP_Hea;
x_7(1) = 0;
T_7(1) = CoolProp.PropsSI('T','P',P_7(1),'Q',x_7(1),fluid);
h_7(1) = CoolProp.PropsSI('H','P',P_7(1),'Q',x_7(1),fluid);
s_7(1) = CoolProp.PropsSI('S','P',P_7(1),'Q',x_7(1),fluid); 

% Energy transfers during the first time step
w_p(1) = h_6(1) - h_5(1);
q_he(1) = h_7(1) - h_6(1);

switch Heat_source
    case 'HPs' % Systems with heat pump to provide the heat at the Heater
        % Hot source temperature (ORES system side of the heat exchanger)
        T_H2(1) = (T_7(1)+T_6(1))/2;
        % T_H2(1) = T(7)
        COP_carnot_HP2(1) = 1/(1-T_amb/T_H2(1));       
        
        if COP_HP2_0 == 0
            COP_HP2(1) = COP_carnot_HP2(1);
        elseif COP_HP2_0 < COP_carnot_HP2(1)
            COP_HP2(1) = COP_HP2_0;
        else
            error('COP_HP2 must be smaller than COP_carnot.')
        end
        w_c_HP2(1) = q_he(1)/COP_HP2(1);
        q_L_HP2(1) = q_he(1) - w_c_HP2(1);
        w_c(1) = w_p(1)+w_c_HP2(1);
%         w_c(1) = w_p(1);
    case 'Free'
        w_c(1) = w_p(1);
end

m_charg(1) = W_c/w_c(1);
% m_charg(1) = W_c/( w_p(1)+w_c_HP2(1) );

dE_p(1) = m_charg(1)*w_p(1)*dt;
dQ_he(1) = m_charg(1)*q_he(1)*dt;

dE_HP2(1) = m_charg(1)*w_c_HP2(1)*dt;

%---------------------------- Dynamic analysis ---------------------------%

Delta_m_ch = 0;
i = 1;
exc_count = 0;
erros = [0, 0];
while Delta_m_ch<Delta_m_disch
    
    i = i+1;    
%     if i == 425
%         i
%     end
    t_charg(i) = t_charg(i-1) + dt;
    Delta_m_ch = Delta_m_ch + m_charg(i-1)*dt;
    m_HPT(n_steps_disch+i) = m_HPT(n_steps_disch+i-1) + m_charg(i-1)*dt;
    rho_HPT(n_steps_disch+i) = m_HPT(n_steps_disch+i)/V_HPT;
    U_HPT(n_steps_disch+i) = U_HPT(n_steps_disch+i-1) + m_charg(i-1)*h_7(i-1)*dt;
    u_HPT(n_steps_disch+i) = U_HPT(n_steps_disch+i)/m_HPT(n_steps_disch+i);

    m_LPT(n_steps_disch+i) = m_LPT(n_steps_disch+i-1) - m_charg(i-1)*dt;
    rho_LPT(n_steps_disch+i) = m_LPT(n_steps_disch+i)/V_LPT;
    U_LPT(n_steps_disch+i) = U_LPT(n_steps_disch+i-1) - m_charg(i-1)*h_5(i-1)*dt;
    u_LPT(n_steps_disch+i) = U_LPT(n_steps_disch+i)/m_LPT(n_steps_disch+i);
    if i==400

    end
    
    try
        if (rho_HPT(n_steps_disch+i)<rho_lim(1) || rho_HPT(n_steps_disch+i)>rho_lim(2)) && (rho_LPT(n_steps_disch+i)<rho_lim(1) || rho_LPT(n_steps_disch+i)>rho_lim(2))
            P_HPT(n_steps_disch+i) = CoolProp.PropsSI('P','U',u_HPT(n_steps_disch+i),'D',rho_HPT(n_steps_disch+i),fluid);
            T_HPT(n_steps_disch+i) = CoolProp.PropsSI('T','U',u_HPT(n_steps_disch+i),'D',rho_HPT(n_steps_disch+i),fluid);
            x_HPT(n_steps_disch+i) = CoolProp.PropsSI('Q','U',u_HPT(n_steps_disch+i),'D',rho_HPT(n_steps_disch+i),fluid);
            h_HPT(n_steps_disch+i) = CoolProp.PropsSI('H','U',u_HPT(n_steps_disch+i),'D',rho_HPT(n_steps_disch+i),fluid);
            s_HPT(n_steps_disch+i) = CoolProp.PropsSI('S','U',u_HPT(n_steps_disch+i),'D',rho_HPT(n_steps_disch+i),fluid);
            ex_HPT(n_steps_disch+i) = (u_HPT(n_steps_disch+i)-u_o) - T_o*(s_HPT(n_steps_disch+i)-s_o) + P_o*(1./rho_HPT(n_steps_disch+i)-1/rho_o);
            EX_HPT(n_steps_disch+i) = m_HPT(n_steps_disch+i)*ex_HPT(n_steps_disch+i);

            P_LPT(n_steps_disch+i) = CoolProp.PropsSI('P','U',u_LPT(n_steps_disch+i),'D',rho_LPT(n_steps_disch+i),fluid);
            T_LPT(n_steps_disch+i) = CoolProp.PropsSI('T','U',u_LPT(n_steps_disch+i),'D',rho_LPT(n_steps_disch+i),fluid);
            x_LPT(n_steps_disch+i) = CoolProp.PropsSI('Q','U',u_LPT(n_steps_disch+i),'D',rho_LPT(n_steps_disch+i),fluid);
            h_LPT(n_steps_disch+i) = CoolProp.PropsSI('H','U',u_LPT(n_steps_disch+i),'D',rho_LPT(n_steps_disch+i),fluid);
            s_LPT(n_steps_disch+i) = CoolProp.PropsSI('S','U',u_LPT(n_steps_disch+i),'D',rho_LPT(n_steps_disch+i),fluid);
            ex_LPT(n_steps_disch+i) = (u_LPT(n_steps_disch+i)-u_o) - T_o*(s_LPT(n_steps_disch+i)-s_o) + P_o*(1./rho_LPT(n_steps_disch+i)-1/rho_o);
            EX_LPT(n_steps_disch+i) = m_LPT(n_steps_disch+i)*ex_LPT(n_steps_disch+i);

        elseif (rho_HPT(n_steps_disch+i)>rho_lim(1) && rho_HPT(n_steps_disch+i)<rho_lim(2)) && (rho_LPT(n_steps_disch+i)<rho_lim(1) || rho_LPT(n_steps_disch+i)>rho_lim(2))
                % If the HPT is inside the unstable region

            if (rho_HPT(n_steps_disch+i)>rho_lim(1) && rho_HPT(n_steps_disch+i)<rho_lim(2)) && (rho_HPT(n_steps_disch+i-1)<rho_lim(1) || rho_HPT(n_steps_disch+i-1)>rho_lim(2))
                % If the HPT was outside the unstable region in the last iteration 
    %             disp('ch. update HPT')

                % Set the inferior limit state
                u_lim_i_HPT = u_HPT(n_steps_disch+i-1);
                rho_lim_i_HPT = rho_HPT(n_steps_disch+i-1);
                P_lim_i_HPT = CoolProp.PropsSI('P','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                T_lim_i_HPT = CoolProp.PropsSI('T','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                x_lim_i_HPT = CoolProp.PropsSI('Q','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                h_lim_i_HPT = CoolProp.PropsSI('H','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                s_lim_i_HPT = CoolProp.PropsSI('S','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);

                % Estimate required mass to exit unstable region
                dm = (rho_lim(2) - rho_HPT(n_steps_disch+i-1))*V_HPT;
                m_lim_s_HPT = m_HPT(n_steps_disch+i-1) + dm;

                % Set the superior limit state
                rho_lim_s_HPT = m_lim_s_HPT/V_HPT;
                U_lim_s_HPT = U_HPT(n_steps_disch+i-1) + dm*h_7(i-1);
                u_lim_s_HPT = U_lim_s_HPT/m_lim_s_HPT;
                P_lim_s_HPT = CoolProp.PropsSI('P','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
                T_lim_s_HPT = CoolProp.PropsSI('T','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
                x_lim_s_HPT = CoolProp.PropsSI('Q','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
                h_lim_s_HPT = CoolProp.PropsSI('H','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
                s_lim_s_HPT = CoolProp.PropsSI('S','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
            end

            % Linearizarion of the unstable region for the HPT
            P_HPT(n_steps_disch+i) = P_lim_i_HPT + (rho_HPT(n_steps_disch+i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(P_lim_s_HPT - P_lim_i_HPT);
            T_HPT(n_steps_disch+i) = T_lim_i_HPT + (rho_HPT(n_steps_disch+i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(T_lim_s_HPT - T_lim_i_HPT);
            x_HPT(n_steps_disch+i) = x_lim_i_HPT + (rho_HPT(n_steps_disch+i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(x_lim_s_HPT - x_lim_i_HPT);
            h_HPT(n_steps_disch+i) = h_lim_i_HPT + (rho_HPT(n_steps_disch+i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(h_lim_s_HPT - h_lim_i_HPT);
            s_HPT(n_steps_disch+i) = s_lim_i_HPT + (rho_HPT(n_steps_disch+i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(s_lim_s_HPT - s_lim_i_HPT);
            ex_HPT(n_steps_disch+i) = (u_HPT(n_steps_disch+i)-u_o) - T_o*(s_HPT(n_steps_disch+i)-s_o) + P_o*(1./rho_HPT(n_steps_disch+i)-1/rho_o);
            EX_HPT(n_steps_disch+i) = m_HPT(n_steps_disch+i)*ex_HPT(n_steps_disch+i);

            P_LPT(n_steps_disch+i) = CoolProp.PropsSI('P','U',u_LPT(n_steps_disch+i),'D',rho_LPT(n_steps_disch+i),fluid);
            T_LPT(n_steps_disch+i) = CoolProp.PropsSI('T','U',u_LPT(n_steps_disch+i),'D',rho_LPT(n_steps_disch+i),fluid);
            x_LPT(n_steps_disch+i) = CoolProp.PropsSI('Q','U',u_LPT(n_steps_disch+i),'D',rho_LPT(n_steps_disch+i),fluid);
            h_LPT(n_steps_disch+i) = CoolProp.PropsSI('H','U',u_LPT(n_steps_disch+i),'D',rho_LPT(n_steps_disch+i),fluid);
            s_LPT(n_steps_disch+i) = CoolProp.PropsSI('S','U',u_LPT(n_steps_disch+i),'D',rho_LPT(n_steps_disch+i),fluid);
            ex_LPT(n_steps_disch+i) = (u_LPT(n_steps_disch+i)-u_o) - T_o*(s_LPT(n_steps_disch+i)-s_o) + P_o*(1./rho_LPT(n_steps_disch+i)-1/rho_o);
            EX_LPT(n_steps_disch+i) = m_LPT(n_steps_disch+i)*ex_LPT(n_steps_disch+i);
        elseif (rho_HPT(n_steps_disch+i)<rho_lim(1) || rho_HPT(n_steps_disch+i)>rho_lim(2)) && (rho_LPT(n_steps_disch+i)>rho_lim(1) && rho_LPT(n_steps_disch+i)<rho_lim(2))
            % If the LPT is inside the unstable region
            P_HPT(n_steps_disch+i) = CoolProp.PropsSI('P','U',u_HPT(n_steps_disch+i),'D',rho_HPT(n_steps_disch+i),fluid);
            T_HPT(n_steps_disch+i) = CoolProp.PropsSI('T','U',u_HPT(n_steps_disch+i),'D',rho_HPT(n_steps_disch+i),fluid);
            x_HPT(n_steps_disch+i) = CoolProp.PropsSI('Q','U',u_HPT(n_steps_disch+i),'D',rho_HPT(n_steps_disch+i),fluid);
            h_HPT(n_steps_disch+i) = CoolProp.PropsSI('H','U',u_HPT(n_steps_disch+i),'D',rho_HPT(n_steps_disch+i),fluid);
            s_HPT(n_steps_disch+i) = CoolProp.PropsSI('S','U',u_HPT(n_steps_disch+i),'D',rho_HPT(n_steps_disch+i),fluid);
            ex_HPT(n_steps_disch+i) = (u_HPT(n_steps_disch+i)-u_o) - T_o*(s_HPT(n_steps_disch+i)-s_o) + P_o*(1./rho_HPT(n_steps_disch+i)-1/rho_o);
            EX_HPT(n_steps_disch+i) = m_HPT(n_steps_disch+i)*ex_HPT(n_steps_disch+i);


            if (rho_LPT(n_steps_disch+i)>rho_lim(1) && rho_LPT(n_steps_disch+i)<rho_lim(2)) && (rho_LPT(n_steps_disch+i-1)<rho_lim(1) || rho_LPT(n_steps_disch+i-1)>rho_lim(2))
                % If the LPT was outside the unstable region in the last iteration 
    %             disp('ch. update LPT')
                % Set the superior limit state
                u_lim_s_LPT = u_LPT(n_steps_disch+i-1);
                rho_lim_s_LPT = rho_LPT(n_steps_disch+i-1);
                P_lim_s_LPT = CoolProp.PropsSI('P','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                T_lim_s_LPT = CoolProp.PropsSI('T','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                x_lim_s_LPT = CoolProp.PropsSI('Q','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                h_lim_s_LPT = CoolProp.PropsSI('H','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                s_lim_s_LPT = CoolProp.PropsSI('S','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);

                % Estimate required mass to exit unstable region
                dm = (rho_LPT(n_steps_disch+i-1) - rho_lim(1))*V_LPT;
                m_lim_i_LPT = m_LPT(n_steps_disch+i-1) - dm;

                % Set the inferior limit state
                rho_lim_i_LPT = m_lim_i_LPT/V_LPT;
                U_lim_i_LPT = U_LPT(n_steps_disch+i-1) - dm*h_5(i-1);
                u_lim_i_LPT = U_lim_i_LPT/m_lim_i_LPT;
                P_lim_i_LPT = CoolProp.PropsSI('P','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
                T_lim_i_LPT = CoolProp.PropsSI('T','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
                x_lim_i_LPT = CoolProp.PropsSI('Q','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
                h_lim_i_LPT = CoolProp.PropsSI('H','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
                s_lim_i_LPT = CoolProp.PropsSI('S','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
            end

            % Linearizarion of the unstable region for the LPT
            P_LPT(n_steps_disch+i) = P_lim_i_LPT + (rho_LPT(n_steps_disch+i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(P_lim_s_LPT - P_lim_i_LPT);
            T_LPT(n_steps_disch+i) = T_lim_i_LPT + (rho_LPT(n_steps_disch+i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(T_lim_s_LPT - T_lim_i_LPT);
            x_LPT(n_steps_disch+i) = x_lim_i_LPT + (rho_LPT(n_steps_disch+i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(x_lim_s_LPT - x_lim_i_LPT);
            h_LPT(n_steps_disch+i) = h_lim_i_LPT + (rho_LPT(n_steps_disch+i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(h_lim_s_LPT - h_lim_i_LPT);
            s_LPT(n_steps_disch+i) = s_lim_i_LPT + (rho_LPT(n_steps_disch+i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(s_lim_s_LPT - s_lim_i_LPT);
            ex_LPT(n_steps_disch+i) = (u_LPT(n_steps_disch+i)-u_o) - T_o*(s_LPT(n_steps_disch+i)-s_o) + P_o*(1./rho_LPT(n_steps_disch+i)-1/rho_o);
            EX_LPT(n_steps_disch+i) = m_LPT(n_steps_disch+i)*ex_LPT(n_steps_disch+i);

        elseif (rho_HPT(n_steps_disch+i)>rho_lim(1) && rho_HPT(n_steps_disch+i)<rho_lim(2)) && (rho_LPT(n_steps_disch+i)>rho_lim(1) && rho_LPT(n_steps_disch+i)<rho_lim(2))
            % If the both tanks are inside the unstable region

            if (rho_HPT(n_steps_disch+i)>rho_lim(1) && rho_HPT(n_steps_disch+i)<rho_lim(2)) && (rho_HPT(n_steps_disch+i-1)<rho_lim(1) || rho_HPT(n_steps_disch+i-1)>rho_lim(2))
                % If the HPT was outside the unstable region in the last iteration 
    %             disp('ch. update HPT 2')
                % Set the inferior limit state
                u_lim_i_HPT = u_HPT(n_steps_disch+i-1);
                rho_lim_i_HPT = rho_HPT(n_steps_disch+i-1);
                P_lim_i_HPT = CoolProp.PropsSI('P','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                T_lim_i_HPT = CoolProp.PropsSI('T','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                x_lim_i_HPT = CoolProp.PropsSI('Q','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                h_lim_i_HPT = CoolProp.PropsSI('H','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);
                s_lim_i_HPT = CoolProp.PropsSI('S','U',u_lim_i_HPT,'D',rho_lim_i_HPT,fluid);

                % Estimate required mass to exit unstable region
                dm = (rho_lim(2) - rho_HPT(n_steps_disch+i-1))*V_HPT;
                m_lim_s_HPT = m_HPT(n_steps_disch+i-1) + dm;

                % Set the superior limit state
                rho_lim_s_HPT = m_lim_s_HPT/V_HPT;
                U_lim_s_HPT = U_HPT(n_steps_disch+i-1) + dm*h_7(i-1);
                u_lim_s_HPT = U_lim_s_HPT/m_lim_s_HPT;
                P_lim_s_HPT = CoolProp.PropsSI('P','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
                T_lim_s_HPT = CoolProp.PropsSI('T','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
                x_lim_s_HPT = CoolProp.PropsSI('Q','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
                h_lim_s_HPT = CoolProp.PropsSI('H','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
                s_lim_s_HPT = CoolProp.PropsSI('S','U',u_lim_s_HPT,'D',rho_lim_s_HPT,fluid);
            end

            % Linearizarion of the unstable region for the HPT
            P_HPT(n_steps_disch+i) = P_lim_i_HPT + (rho_HPT(n_steps_disch+i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(P_lim_s_HPT - P_lim_i_HPT);
            T_HPT(n_steps_disch+i) = T_lim_i_HPT + (rho_HPT(n_steps_disch+i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(T_lim_s_HPT - T_lim_i_HPT);
            x_HPT(n_steps_disch+i) = x_lim_i_HPT + (rho_HPT(n_steps_disch+i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(x_lim_s_HPT - x_lim_i_HPT);
            h_HPT(n_steps_disch+i) = h_lim_i_HPT + (rho_HPT(n_steps_disch+i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(h_lim_s_HPT - h_lim_i_HPT);
            s_HPT(n_steps_disch+i) = s_lim_i_HPT + (rho_HPT(n_steps_disch+i)-rho_lim_i_HPT)/(rho_lim_s_HPT-rho_lim_i_HPT)*(s_lim_s_HPT - s_lim_i_HPT);
            ex_HPT(n_steps_disch+i) = (u_HPT(n_steps_disch+i)-u_o) - T_o*(s_HPT(n_steps_disch+i)-s_o) + P_o*(1./rho_HPT(n_steps_disch+i)-1/rho_o);
            EX_HPT(n_steps_disch+i) = m_HPT(n_steps_disch+i)*ex_HPT(n_steps_disch+i);

            if (rho_LPT(n_steps_disch+i)>rho_lim(1) && rho_LPT(n_steps_disch+i)<rho_lim(2)) && (rho_LPT(n_steps_disch+i-1)<rho_lim(1) || rho_LPT(n_steps_disch+i-1)>rho_lim(2))
                % If the LPT was outside the unstable region in the last iteration 
    %             disp('ch. update LPT 2')
                % Set the superior limit state
                u_lim_s_LPT = u_LPT(n_steps_disch+i-1);
                rho_lim_s_LPT = rho_LPT(n_steps_disch+i-1);
                P_lim_s_LPT = CoolProp.PropsSI('P','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                T_lim_s_LPT = CoolProp.PropsSI('T','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                x_lim_s_LPT = CoolProp.PropsSI('Q','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                h_lim_s_LPT = CoolProp.PropsSI('H','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);
                s_lim_s_LPT = CoolProp.PropsSI('S','U',u_lim_s_LPT,'D',rho_lim_s_LPT,fluid);

                % Estimate required mass to exit unstable region
                dm = (rho_LPT(n_steps_disch+i-1) - rho_lim(1))*V_LPT;
                m_lim_i_LPT = m_LPT(n_steps_disch+i-1) - dm;

                % Set the inferior limit state
                rho_lim_i_LPT = m_lim_i_LPT/V_LPT;
                U_lim_i_LPT = U_LPT(n_steps_disch+i-1) - dm*h_5(i-1);
                u_lim_i_LPT = U_lim_i_LPT/m_lim_i_LPT;
                P_lim_i_LPT = CoolProp.PropsSI('P','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
                T_lim_i_LPT = CoolProp.PropsSI('T','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
                x_lim_i_LPT = CoolProp.PropsSI('Q','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
                h_lim_i_LPT = CoolProp.PropsSI('H','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
                s_lim_i_LPT = CoolProp.PropsSI('S','U',u_lim_i_LPT,'D',rho_lim_i_LPT,fluid);
            end

            % Linearizarion of the unstable region for the LPT
            P_LPT(n_steps_disch+i) = P_lim_i_LPT + (rho_LPT(n_steps_disch+i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(P_lim_s_LPT - P_lim_i_LPT);
            T_LPT(n_steps_disch+i) = T_lim_i_LPT + (rho_LPT(n_steps_disch+i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(T_lim_s_LPT - T_lim_i_LPT);
            x_LPT(n_steps_disch+i) = x_lim_i_LPT + (rho_LPT(n_steps_disch+i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(x_lim_s_LPT - x_lim_i_LPT);
            h_LPT(n_steps_disch+i) = h_lim_i_LPT + (rho_LPT(n_steps_disch+i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(h_lim_s_LPT - h_lim_i_LPT);
            s_LPT(n_steps_disch+i) = s_lim_i_LPT + (rho_LPT(n_steps_disch+i)-rho_lim_i_LPT)/(rho_lim_s_LPT-rho_lim_i_LPT)*(s_lim_s_LPT - s_lim_i_LPT);
            ex_LPT(n_steps_disch+i) = (u_LPT(n_steps_disch+i)-u_o) - T_o*(s_LPT(n_steps_disch+i)-s_o) + P_o*(1./rho_LPT(n_steps_disch+i)-1/rho_o);
            EX_LPT(n_steps_disch+i) = m_LPT(n_steps_disch+i)*ex_LPT(n_steps_disch+i);
        end
        exc_count = 0;
    catch ME
%         ME
%         getReport(ME)
        erros(i,:) = [1,i];             
               
%         ME.identifier
        switch ME.identifier
            case 'SWIG:RuntimeError'
                P_HPT(n_steps_disch+i) = P_HPT(n_steps_disch+i-1) + (rho_HPT(n_steps_disch+i)-rho_HPT(n_steps_disch+i-1))/(rho_HPT(n_steps_disch+i-1)-rho_HPT(n_steps_disch+i-2))*(P_HPT(n_steps_disch+i-1) - P_HPT(n_steps_disch+i-2));
                T_HPT(n_steps_disch+i) = T_HPT(n_steps_disch+i-1) + (rho_HPT(n_steps_disch+i)-rho_HPT(n_steps_disch+i-1))/(rho_HPT(n_steps_disch+i-1)-rho_HPT(n_steps_disch+i-2))*(T_HPT(n_steps_disch+i-1) - T_HPT(n_steps_disch+i-2));
                x_HPT(n_steps_disch+i) = x_HPT(n_steps_disch+i-1) + (rho_HPT(n_steps_disch+i)-rho_HPT(n_steps_disch+i-1))/(rho_HPT(n_steps_disch+i-1)-rho_HPT(n_steps_disch+i-2))*(x_HPT(n_steps_disch+i-1) - x_HPT(n_steps_disch+i-2));
                h_HPT(n_steps_disch+i) = h_HPT(n_steps_disch+i-1) + (rho_HPT(n_steps_disch+i)-rho_HPT(n_steps_disch+i-1))/(rho_HPT(n_steps_disch+i-1)-rho_HPT(n_steps_disch+i-2))*(h_HPT(n_steps_disch+i-1) - h_HPT(n_steps_disch+i-2));
                s_HPT(n_steps_disch+i) = s_HPT(n_steps_disch+i-1) + (rho_HPT(n_steps_disch+i)-rho_HPT(n_steps_disch+i-1))/(rho_HPT(n_steps_disch+i-1)-rho_HPT(n_steps_disch+i-2))*(s_HPT(n_steps_disch+i-1) - s_HPT(n_steps_disch+i-2));
                
                ex_HPT(n_steps_disch+i) = (u_HPT(n_steps_disch+i)-u_o) - T_o*(s_HPT(n_steps_disch+i)-s_o) + P_o*(1./rho_HPT(n_steps_disch+i)-1/rho_o);
                EX_HPT(n_steps_disch+i) = m_HPT(n_steps_disch+i)*ex_HPT(n_steps_disch+i);
                
                P_LPT(n_steps_disch+i) = P_LPT(n_steps_disch+i-1) + (rho_LPT(n_steps_disch+i)-rho_LPT(n_steps_disch+i-1))/(rho_LPT(n_steps_disch+i-1)-rho_LPT(n_steps_disch+i-2))*(P_LPT(n_steps_disch+i-1) - P_LPT(n_steps_disch+i-2));
                T_LPT(n_steps_disch+i) = T_LPT(n_steps_disch+i-1) + (rho_LPT(n_steps_disch+i)-rho_LPT(n_steps_disch+i-1))/(rho_LPT(n_steps_disch+i-1)-rho_LPT(n_steps_disch+i-2))*(T_LPT(n_steps_disch+i-1) - T_LPT(n_steps_disch+i-2));
                x_LPT(n_steps_disch+i) = x_LPT(n_steps_disch+i-1) + (rho_LPT(n_steps_disch+i)-rho_LPT(n_steps_disch+i-1))/(rho_LPT(n_steps_disch+i-1)-rho_LPT(n_steps_disch+i-2))*(x_LPT(n_steps_disch+i-1) - x_LPT(n_steps_disch+i-2));
                h_LPT(n_steps_disch+i) = h_LPT(n_steps_disch+i-1) + (rho_LPT(n_steps_disch+i)-rho_LPT(n_steps_disch+i-1))/(rho_LPT(n_steps_disch+i-1)-rho_LPT(n_steps_disch+i-2))*(h_LPT(n_steps_disch+i-1) - h_LPT(n_steps_disch+i-2));
                s_LPT(n_steps_disch+i) = s_LPT(n_steps_disch+i-1) + (rho_LPT(n_steps_disch+i)-rho_LPT(n_steps_disch+i-1))/(rho_LPT(n_steps_disch+i-1)-rho_LPT(n_steps_disch+i-2))*(s_LPT(n_steps_disch+i-1) - s_LPT(n_steps_disch+i-2));
                
                ex_LPT(n_steps_disch+i) = (u_LPT(n_steps_disch+i)-u_o) - T_o*(s_LPT(n_steps_disch+i)-s_o) + P_o*(1./rho_LPT(n_steps_disch+i)-1/rho_o);
                EX_LPT(n_steps_disch+i) = m_LPT(n_steps_disch+i)*ex_LPT(n_steps_disch+i);
                
                exc_count = exc_count + 1;
                if exc_count >= 5
                   error('5 consecutive errors.') 
                end
            otherwise
                rethrow(ME)
        end
%         i
        
%         if t_generation < 3600
%             break;
%         end
    end
        % Cycle state update
    % State at LPT outlet / Pump inlet
        P_5(i) = P_LPT(n_steps_disch+i); % INCLUIR VERIFICAÇÃO DE ESTADO LIQUIDO SATURADO OU SUBRESFRIADO
        x_5(i) = 0;
        T_5(i) = CoolProp.PropsSI('T','P',P_5(i),'Q',x_5(i),fluid);
        h_5(i) = CoolProp.PropsSI('H','P',P_5(i),'Q',x_5(i),fluid);
        s_5(i) = CoolProp.PropsSI('S','P',P_5(i),'Q',x_5(i),fluid);

    % State at Pump outlet / Heater inlet
    %     P_6(i) = P_HPT(1); % VERIFICAR
        P_6(i) = P_HPT(n_steps_disch+i);

        h_p(i) = CoolProp.PropsSI('H','P',P_6(i),'S',s_5(i),fluid);
        h_6(i) = pump('h_in',h_5(i),'h_s',h_p(i),'eta_p',eta_p);
        T_6(i) = CoolProp.PropsSI('T','P',P_6(i),'H',h_6(i),fluid);
        s_6(i) = CoolProp.PropsSI('S','P',P_6(i),'H',h_6(i),fluid);    

    % State at Heater outlet / HPT inlet
        P_7(i) = P_6(i) - DeltaP_Hea;
        x_7(i) = 0;
        T_7(i) = CoolProp.PropsSI('T','P',P_7(i),'Q',x_7(i),fluid);
        h_7(i) = CoolProp.PropsSI('H','P',P_7(i),'Q',x_7(i),fluid);
        s_7(i) = CoolProp.PropsSI('S','P',P_7(i),'Q',x_7(i),fluid);

    % Energy transfers during i'th time step
        w_p(i) = h_6(i) - h_5(i);
        q_he(i) = h_7(i) - h_6(i);

        switch Heat_source
            case 'HPs' % Systems with heat pump to provide the heat at the Heater
                % Hot source temperature (ORES system side of the heat exchanger)
                T_H2(i) = (T_7(i)+T_6(i))/2;
                % T_H2(1) = T(7)
                COP_carnot_HP2(i) = 1/(1-T_amb/T_H2(i));       

                if COP_HP2_0 == 0
                    COP_HP2(i) = COP_carnot_HP2(i);
                elseif COP_HP2_0 < COP_carnot_HP2(i)
                    COP_HP2(i) = COP_HP2_0;
                else
                    error('COP_HP2 must be smaller than COP_carnot.')
                end
                w_c_HP2(i) = q_he(i)/COP_HP2(i);
                q_L_HP2(i) = q_he(i) - w_c_HP2(i);
                w_c(i) = w_p(i)+w_c_HP2(i);
    %             w_c(i) = w_p(i);
            case 'Free'
                w_c(i) = w_p(i);
        end

        m_charg(i) = W_c/w_c(i);
        dE_p(i) = m_charg(i)*w_p(i)*dt;
        dQ_he(i) = m_charg(i)*q_he(i)*dt;

        dE_HP2(i) = m_charg(i)*w_c_HP2(i)*dt;

        if x_HPT(n_steps_disch+i)<x_HPT_f_ch
            disp('x_HPT < 0.02 during charge.')
            break;
        elseif x_LPT(n_steps_disch+i)>=x_LPT_f_ch
            disp('x_LPT > 0.98 during charge.')
            break;
        end
    
end
  
%---------------------------- Results display ----------------------------%
if plotBool
    % Properties over time
    figure('color',[1 1 1])
    subplot(3,2,1)
    plot(t_charg,m_charg,'k')
    ylabel('mass flow rate [kg/s]')
    xtickformat('%,0.0f')
    ytickformat('%0.1f')
    grid on;
    
    subplot(3,2,2) % Total mass in the tanks 
    plot(t_charg,m_HPT(n_steps_disch+1:end),'k'); hold on;
    plot(t_charg,m_LPT(n_steps_disch+1:end),'color',[0.7 0.7 0.7]); 
    grid on;
    ylabel('Total mass [kg]')
    xtickformat('%,0.0f')
    ytickformat('%0.0f')
    legend('m_{HPT}','m_{LPT}','location','east')
    
    subplot(3,2,3) % Pressure at the tanks
    plot(t_charg,P_HPT(n_steps_disch+1:end)./1000,'k'); hold on;
    plot(t_charg,P_LPT(n_steps_disch+1:end)./1000,'color',[0.7 0.7 0.7]); 
    grid on;
    legend('P_{HPT}','P_{LPT}','location','east')
    ylabel('Pressure [kPa]')
    xtickformat('%,0.0f')
    ytickformat('%,0.0f')
    
    subplot(3,2,4) % Temperature at the tanks
    plot(t_charg,T_HPT(n_steps_disch+1:end),'k'); hold on;
    plot(t_charg,T_LPT(n_steps_disch+1:end),'color',[0.7 0.7 0.7]); 
    grid on;
    ylabel('Temperature [K]')
    legend('T_{HPT}','T_{LPT}','location','east')
    xtickformat('%,0.0f')
    
    subplot(3,2,5) % Quality of the fluid in the tanks
    plot(t_charg,x_HPT(n_steps_disch+1:end),'k'); hold on;
    plot(t_charg,x_LPT(n_steps_disch+1:end),'color',[0.7 0.7 0.7]); 
    grid on;
    ylim([0 1])
    ylabel('Quality [-]')
    legend('x_{HPT}','x_{LPT}','location','northwest')
    xtickformat('%,0.0f')
    ytickformat('%0.1f')
    
    subplot(3,2,6)
    plot(t_charg,rho_HPT(n_steps_disch+1:end),'k'); hold on; 
    plot(t_charg,rho_LPT(n_steps_disch+1:end),'color',[0.7 0.7 0.7]); 
    grid on;
    ylabel('Density [kg/m^3]')
    xlabel('Time [s]')
    legend('\rho_{HPT}','\rho_{LPT}','location','east')
    xtickformat('%,0.0f')
    ytickformat('%0.0f')
    
    % Plot the T-s diagramm for the fluid with the defined isobaric lines
%     TSDiag=plotTS(fluid,'Iso_P',[P_HPT(n+1),P_LPT(n+1),P_HPT(end),P_LPT(end)]);
    plotTS(fluid,'Iso_P',[P_HPT(n_steps_disch+1),P_LPT(n_steps_disch+1),P_HPT(end),P_LPT(end)]);
    plot(s_HPT(n_steps_disch+1:end)./1000,T_HPT(n_steps_disch+1:end),'k','LineWidth',2)
    handle_charg(1) = plot(s_LPT(n_steps_disch+1:end)./1000,T_LPT(n_steps_disch+1:end),'k','LineWidth',2);
    
    plotTS(fluid,'Process','Pump_comp',P_5(1),h_5(1),P_6(1),eta_p,'plotColor','b');%'k');
    handle_charg(2) = plotTS(fluid,'Process','Isobaric',P_6(1),s_6(1),s_7(1),'plotColor','b');%'k');
    
    plotTS(fluid,'Process','Pump_comp',P_5(end),h_5(end),P_6(end),eta_p,'plotColor','r')%[0.7 0.7 0.7]);
    handle_charg(3) = plotTS(fluid,'Process','Isobaric',P_6(end),s_6(end),s_7(end),'plotColor','r');%[0.7 0.7 0.7]);
    grid on;
    % T-s diagramm during charging phase
    ylim([T_LPT(1)-10 T_2(1)+20]);
    text([s_5(1);s_6(1)-15;s_7(1)]./1000,[T_5(1)-5;T_6(1)+10;T_7(end)],{'5','6','7'},'FontSize',12)
%     text([s_5(end);s_6(end);s_7(end)]./1000,...
%         [T_5(end);T_6(end);T_7(end)],{'5\prime','6\prime','7\prime'},...
%         'Color',[0.5 0.5 0.5]);
    plot([s_5(1);s_6(1);s_7(1)]./1000,[T_5(1);T_6(1);T_7(1)],'b.','MarkerSize',20)
    plot([s_5(end);s_6(end);s_7(end)]./1000,...
         [T_5(end);T_6(end);T_7(end)],'k.','MarkerSize',20,'Color','r')
    legend(handle_charg,{'State at the tanks','Process diagramm at t=0 s',...
        'Process diagramm at t=t_{charge}'},'location','northwest')
    xtickformat('%0.1f')
    fig = gcf;
    fig.PaperOrientation = 'landscape';
    fig.PaperSize = [6 4.5];
    %xlim([floor(s_LPT(1)./1000) ceil(s_3(end)./100)./10])
end

%=================== Full-cycle performance parameters ===================%
    E_p = sum(dE_p); % Total work consumed in the pump, in J
    E_t = sum(dE_t); % Total work generated in the turbine, in J
%     Q_in = m_charg*sum(q_he*dt)+sum(m_disch.*q_ev*dt); % Total consumed heat

    E_HP2 = sum(dE_HP2);
    
    E_c = E_p + E_HP2;
    E_g = sum(dE_net);
    eta_RT = E_g/E_c;
    % Key performance indexes
%     eta_RT = E_g/E_c;
%     eta_I = (E_g-E_c)/Q_in;
%     eta_II = E_g/(Q_in*(1-T_amb/T_H)+E_c);
    rho_E_L = ( E_t/3600 )/( 1000*(V_HPT + V_LPT) ); % Energy density [Wh/L]
    %rho_E_kg = ; % Energy density [Wh/kg]
    
    if sum(erros(:,1))>0
        erro = 1;
%         eta_RT = NaN;
%         rho_E_L = NaN;
    end
    
    % Exergy at the tanks over time
    v_HPT = 1./rho_HPT;
    v_LPT = 1./rho_LPT;    
    ex_HPT = u_HPT - u_o + P_o.*(1./rho_HPT-v_o) - T_o.*(s_HPT - s_o);
    ex_LPT = u_LPT - u_o + P_o.*(1./rho_LPT-v_o) - T_o.*(s_LPT - s_o);
    
    Ex_HPT_kJ_kg = ex_HPT./1000;
    Ex_HPT_kJ_l = ((ex_HPT)./1000)./(v_HPT*1000);                % [kJ/L]
    Ex_HPT_Wh_kg = ex_HPT./3600;
    Ex_HPT_Wh_l = (ex_HPT./3600)./(v_HPT*1000);
    
    Ex_LPT_kJ_kg = ex_LPT./1000;
    Ex_LPT_kJ_l = ((ex_LPT)./1000)./(v_LPT*1000);                % [kJ/L]
    Ex_LPT_Wh_kg = ex_LPT./3600;
    Ex_LPT_Wh_l = (ex_LPT./3600)./(v_LPT*1000);
    
    % CAPEX    
    C_fluid = fluid_cost(fluid,m_HPT(1) + m_LPT(1));
    
%     C_tank_H = vessel_cost('Bare Module',max(P_HPT),V_HPT,'Vertical','Carbon Steel')
%     C_tank_L = vessel_cost('Bare Module',max(P_LPT),V_LPT,'Vertical','Carbon Steel')
%     C_tank = C_tank_H + C_tank_L;
    
%     Storage_Cost_per_m3 = 750; % Pressure tank cost [$/m3]
%     C_tank = Storage_Cost_per_m3*(V_HPT + V_LPT) ;
    C_tank = vessel_cost('Matheus',1,(V_HPT + V_LPT),'','');

    C_turb = turbine_cost('Bare Module',max(m_disch.*w_t)/1000,'Axial Gas Turbine','Stainless Steel');
%     C_turb = turbine_cost('six-tenth',max(m_disch.*w_t)/1000,'Axial Gas Turbine','Carbon Steel');
    
    C_pump = pump_cost(max(m_charg.*w_p)/1000,P_HPT(1),'Centrifugal','Carbon Steel');
    
    C_aux = cost_index(2018)/cost_index(2009)*328*( (sum(m_disch.*w_t.*dt)./3600) )./1000;
    
    CAPEX = C_fluid + C_tank + C_turb + C_pump + C_aux;
    [V_HPT;V_LPT];
end