function [eta_RT,eta_I,eta_II,q_in,w_net]=ORES_st(varargin)
% ORES_ST ORES steady-state operation model
%   ORES_st()
% sets of inputs and outputs:                                             %
%                                                                         %
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
% P_amb = 101325;

% System parameters and requirements
fluid = 'R141b';
% W_t = 2000;
% W_p = 0.2*W_t;

% Assumptions
eta_t = 0.80;
eta_p = 0.75;

P_HPT_default = 3500000;
x_HPT_default = 0.02;

% P_LPT_default = 1000000;
T_LPT_default = T_amb;
x_LPT_default = 0.98;
DT_SH = 0.01;

Heat_source = 'HPs';
COP_HP1_0 = 0;
COP_HP2_0 = 0;
%------- Basic parameters: Fluid definition and ambient conditions -------%
    if any(strcmp('fluid',varargin))
        fluid = varargin{find(strcmp('fluid',varargin))+1};
        if ~ischar(fluid)
            error('fluid must be a string, not a %s.',class(fluid))
        end
    end
    if any(strcmp('T_amb',varargin))
        T_amb = varargin{find(strcmp('T_amb',varargin))+1};
        if ~isnumeric(T_amb)
            error('T_amb must be a number, not a %s.',class(T_amb))
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
    
%------------------------------ State at HPT -----------------------------%
    if any(strcmp('P_HPT',varargin))    % INPUT P_HPT
        if ~isnumeric(varargin{find(strcmp('P_HPT',varargin))+1})
            error('P_HPT must be a number, not a %s.',class(varargin{find(strcmp('P_HPT',varargin))+1}))
        end
        P_HPT = varargin{find(strcmp('P_HPT',varargin))+1};
        if any(strcmp('x_HPT',varargin))    % INPUT P_HPT + x_HPT
            if ~isnumeric(varargin{find(strcmp('x_HPT',varargin))+1})
                error('x_HPT must be a number, not a %s.',class(varargin{find(strcmp('x_HPT',varargin))+1}))
            elseif any(strcmp('T_HPT',varargin))
                warning('P_HPT and x_HPT are already determined, given T_HPT value will be ignored.')
            end
            x_HPT = varargin{find(strcmp('x_HPT',varargin))+1};
            T_HPT = CoolProp.PropsSI('T','P',P_HPT,'Q',x_HPT,fluid);
            u_HPT = CoolProp.PropsSI('U','P',P_HPT,'Q',x_HPT,fluid);
            h_HPT = CoolProp.PropsSI('H','P',P_HPT,'Q',x_HPT,fluid);
            s_HPT = CoolProp.PropsSI('S','P',P_HPT,'Q',x_HPT,fluid);
            ph_HPT = CoolProp.PropsSI('Phase','P',P_HPT,'Q',x_HPT,fluid);
        elseif any(strcmp('T_HPT',varargin))    % INPUT P_HPT + T_HPT
            if ~isnumeric(varargin{find(strcmp('T_HPT',varargin))+1})
                error('T_HPT must be a number, not a %s.',class(varargin{find(strcmp('T_HPT',varargin))+1}))
            end
            T_HPT = varargin{find(strcmp('T_HPT',varargin))+1};
            %x_HPT = CoolProp.PropsSI('Q','P',P_HPT,'T',T_HPT,fluid);
            u_HPT = CoolProp.PropsSI('U','P',P_HPT,'T',T_HPT,fluid);
            h_HPT = CoolProp.PropsSI('H','P',P_HPT,'T',T_HPT,fluid);
            s_HPT = CoolProp.PropsSI('S','P',P_HPT,'T',T_HPT,fluid);
            ph_HPT = CoolProp.PropsSI('Phase','P',P_HPT,'T',T_HPT,fluid);
        else                % INPUT P_HPT + DEFAULT
            x_HPT = x_HPT_default;
            T_HPT = CoolProp.PropsSI('T','P',P_HPT,'Q',x_HPT,fluid);
            u_HPT = CoolProp.PropsSI('U','P',P_HPT,'Q',x_HPT,fluid);
            h_HPT = CoolProp.PropsSI('H','P',P_HPT,'Q',x_HPT,fluid);
            s_HPT = CoolProp.PropsSI('S','P',P_HPT,'Q',x_HPT,fluid);
            ph_HPT = CoolProp.PropsSI('Phase','P',P_HPT,'Q',x_HPT,fluid);
        end
    else
        if any(strcmp('x_HPT',varargin))    % INPUT x_HPT
            if ~isnumeric(varargin{find(strcmp('x_HPT',varargin))+1})
                error('x_HPT must be a number, not a %s.',class(varargin{find(strcmp('x_HPT',varargin))+1}))
            end
            x_HPT = varargin{find(strcmp('x_HPT',varargin))+1};
            if any(strcmp('T_HPT',varargin))  % INPUT x_HPT + T_HPT
                if ~isnumeric(varargin{find(strcmp('T_HPT',varargin))+1})
                    error('T_HPT must be a number, not a %s.',class(varargin{find(strcmp('T_HPT',varargin))+1}))
                end
                T_HPT = varargin{find(strcmp('T_HPT',varargin))+1};
                P_HPT = CoolProp.PropsSI('P','Q',x_HPT,'T',T_HPT,fluid);
                rho_HPT = CoolProp.PropsSI('D','Q',x_HPT,'T',T_HPT,fluid);
                u_HPT = CoolProp.PropsSI('U','Q',x_HPT,'T',T_HPT,fluid);
                h_HPT = CoolProp.PropsSI('H','Q',x_HPT,'T',T_HPT,fluid);
                s_HPT = CoolProp.PropsSI('S','Q',x_HPT,'T',T_HPT,fluid);
                ph_HPT = CoolProp.PropsSI('Phase','T',T_HPT,'Q',x_HPT,fluid);
            else                            % INPUT x_HPT + DEFAULT
                P_HPT = P_HPT_default;
                T_HPT = CoolProp.PropsSI('T','Q',x_HPT,'P',P_HPT,fluid);
                rho_HPT = CoolProp.PropsSI('D','Q',x_HPT,'P',P_HPT,fluid);
                u_HPT = CoolProp.PropsSI('U','Q',x_HPT,'P',P_HPT,fluid);
                h_HPT = CoolProp.PropsSI('H','Q',x_HPT,'P',P_HPT,fluid);
                s_HPT = CoolProp.PropsSI('S','Q',x_HPT,'P',P_HPT,fluid);
                ph_HPT = CoolProp.PropsSI('Phase','Q',x_HPT,'T',T_HPT,fluid);
            end
        elseif any(strcmp('T_HPT',varargin))  % INPUT T_HPT + DEFAULT
            if ~isnumeric(varargin{find(strcmp('T_HPT',varargin))+1})
                error('T_HPT must be a number, not a %s.',class(varargin{find(strcmp('T_HPT',varargin))+1}))
            end
            T_HPT = varargin{find(strcmp('T_HPT',varargin))+1};
            x_HPT = x_HPT_default;  
            P_HPT = CoolProp.PropsSI('P','Q',x_HPT,'T',T_HPT,fluid);
            rho_HPT = CoolProp.PropsSI('D','Q',x_HPT,'T',T_HPT,fluid);
            u_HPT = CoolProp.PropsSI('U','Q',x_HPT,'T',T_HPT,fluid);
            h_HPT = CoolProp.PropsSI('H','Q',x_HPT,'T',T_HPT,fluid);
            s_HPT = CoolProp.PropsSI('S','Q',x_HPT,'T',T_HPT,fluid);
            ph_HPT = CoolProp.PropsSI('Phase','P',P_HPT,'Q',x_HPT,fluid);
        else                                % If no variable for HPT is available
            warning('No input for HPT, default values of x=0 and P=2.131 kPa assumed.')
            P_HPT = P_HPT_default;
            x_HPT = x_HPT_default;
            T_HPT = CoolProp.PropsSI('T','Q',x_HPT,'P',P_HPT,fluid);
            rho_HPT = CoolProp.PropsSI('D','Q',x_HPT,'P',P_HPT,fluid);
            u_HPT = CoolProp.PropsSI('U','Q',x_HPT,'P',P_HPT,fluid);
            h_HPT = CoolProp.PropsSI('H','Q',x_HPT,'P',P_HPT,fluid);
            s_HPT = CoolProp.PropsSI('S','Q',x_HPT,'P',P_HPT,fluid);
            ph_HPT = CoolProp.PropsSI('Phase','Q',x_HPT,'P',P_HPT,fluid);
        end
    end

    
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
        T_H = varargin{find(strcmp('T_H',varargin))+1};
        if ~isnumeric(varargin{find(strcmp('T_2',varargin))+1})
            error('T_2 must be a number, not a %s.',class(varargin{find(strcmp('T_2',varargin))+1}))
        elseif ~isnumeric(varargin{find(strcmp('T_H',varargin))+1})
            error('T_H must be a number, not a %s.',class(varargin{find(strcmp('T_H',varargin))+1}))
        end
    elseif any(strcmp('T_2',varargin))
        if ~isnumeric(varargin{find(strcmp('T_2',varargin))+1})
            error('T_2 must be a number, not a %s.',class(varargin{find(strcmp('T_2',varargin))+1}))
        end
        T_2 = varargin{find(strcmp('T_2',varargin))+1};
        T_H = T_2 + 20;
    elseif any(strcmp('T_H',varargin))
        if ~isnumeric(varargin{find(strcmp('T_H',varargin))+1})
            error('T_H must be a number, not a %s.',class(varargin{find(strcmp('T_H',varargin))+1}))
        end
        T_H = varargin{find(strcmp('T_H',varargin))+1};
        T_2 = T_H-20;
    elseif any(strcmp('DT_SH',varargin))
        if ~isnumeric(varargin{find(strcmp('DT_SH',varargin))+1})
            error('DT_SH must be a number, not a %s.',class(varargin{find(strcmp('DT_SH',varargin))+1}))
        end
        DT_SH = varargin{find(strcmp('DT_SH',varargin))+1};
        T_sat_Ev = CoolProp.PropsSI('T','P',P_HPT,'Q',0,fluid);
        T_2 = T_sat_Ev + DT_SH;
        T_H = T_2 + 20;
    else
        T_sat_Ev = CoolProp.PropsSI('T','P',P_HPT,'Q',0,fluid);
        T_2 = T_sat_Ev + DT_SH;
        T_H = T_2 + 20;
    end            
    
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
    else
        DeltaP_Con = 0;
    end
    
%----------------------------- State at LPT ------------------------------%
    if any(strcmp('P_LPT',varargin))
        if ~isnumeric(varargin{find(strcmp('P_LPT',varargin))+1})
            error('P_LPT must be a number, not a %s.',class(varargin{find(strcmp('P_LPT',varargin))+1}))
        end
        P_LPT = varargin{find(strcmp('P_LPT',varargin))+1};
        if any(strcmp('x_LPT',varargin))      % If P and x are available
            if ~isnumeric(varargin{find(strcmp('x_LPT',varargin))+1})
                error('x_LPT must be a number, not a %s.',class(varargin{find(strcmp('x_LPT',varargin))+1}))
            elseif any(strcmp('T_LPT',varargin))
                warning('P_LPT and x_LPT are already determined, given T_LPT value will be ignored.')
            end
            x_LPT = varargin{find(strcmp('x_LPT',varargin))+1};
            T_LPT = CoolProp.PropsSI('T','P',P_LPT,'Q',x_LPT,fluid);
            rho_LPT = CoolProp.PropsSI('D','P',P_LPT,'Q',x_LPT,fluid);
            u_LPT = CoolProp.PropsSI('U','P',P_LPT,'Q',x_LPT,fluid);
            h_LPT = CoolProp.PropsSI('H','P',P_LPT,'Q',x_LPT,fluid);
            s_LPT = CoolProp.PropsSI('S','P',P_LPT,'Q',x_LPT,fluid);
            ph_LPT = CoolProp.PropsSI('Phase','P',P_LPT,'Q',x_LPT,fluid);
        elseif any(strcmp('T_LPT',varargin))  % If P and x are available
            if ~isnumeric(varargin{find(strcmp('T_LPT',varargin))+1})
                error('T_LPT must be a number, not a %s.',class(varargin{find(strcmp('T_LPT',varargin))+1}))
            end
            T_LPT = varargin{find(strcmp('T_LPT',varargin))+1};
            x_LPT = CoolProp.PropsSI('Q','P',P_LPT,'T',T_LPT,fluid);
            rho_LPT = CoolProp.PropsSI('D','P',P_LPT,'T',T_LPT,fluid);
            u_LPT = CoolProp.PropsSI('U','P',P_LPT,'T',T_LPT,fluid);
            h_LPT = CoolProp.PropsSI('H','P',P_LPT,'T',T_LPT,fluid);            
            s_LPT = CoolProp.PropsSI('S','P',P_LPT,'T',T_LPT,fluid);
            ph_LPT = CoolProp.PropsSI('Phase','P',P_LPT,'T',T_LPT,fluid);
        else                                % If only P is available
            x_LPT = x_LPT_default;
            T_LPT = CoolProp.PropsSI('T','P',P_LPT,'Q',x_LPT,fluid);
            rho_LPT = CoolProp.PropsSI('D','P',P_LPT,'Q',x_LPT,fluid);
            u_LPT = CoolProp.PropsSI('U','P',P_LPT,'Q',x_LPT,fluid);
            h_LPT = CoolProp.PropsSI('H','P',P_LPT,'Q',x_LPT,fluid);
            s_LPT = CoolProp.PropsSI('S','P',P_LPT,'Q',x_LPT,fluid);
            ph_LPT = CoolProp.PropsSI('Phase','P',P_LPT,'Q',x_LPT,fluid);
        end 
    else
        if any(strcmp('x_LPT',varargin))
            if ~isnumeric(varargin{find(strcmp('x_LPT',varargin))+1})
                error('x_LPT must be a number, not a %s.',class(varargin{find(strcmp('x_LPT',varargin))+1}))
            end
            x_LPT = varargin{find(strcmp('x_LPT',varargin))+1};
            if any(strcmp('T_LPT',varargin))  % If x and T are available
                if ~isnumeric(varargin{find(strcmp('T_LPT',varargin))+1})
                    error('T_LPT must be a number, not a %s.',class(varargin{find(strcmp('T_LPT',varargin))+1}))
                end
                T_LPT = varargin{find(strcmp('T_LPT',varargin))+1};
                P_LPT = CoolProp.PropsSI('P','T',T_LPT,'Q',x_LPT,fluid);
                rho_LPT = CoolProp.PropsSI('D','T',T_LPT,'Q',x_LPT,fluid);
                u_LPT = CoolProp.PropsSI('U','T',T_LPT,'Q',x_LPT,fluid);
                h_LPT = CoolProp.PropsSI('H','T',T_LPT,'Q',x_LPT,fluid);
                s_LPT = CoolProp.PropsSI('S','T',T_LPT,'Q',x_LPT,fluid);
                ph_LPT = CoolProp.PropsSI('Phase','T',T_LPT,'Q',x_LPT,fluid);
            else                            % If only x is available
                T_LPT = T_LPT_default;
                P_LPT = CoolProp.PropsSI('P','T',T_LPT,'Q',x_LPT,fluid);
                rho_LPT = CoolProp.PropsSI('D','T',T_LPT,'Q',x_LPT,fluid);
                u_LPT = CoolProp.PropsSI('U','T',T_LPT,'Q',x_LPT,fluid);
                h_LPT = CoolProp.PropsSI('H','T',T_LPT,'Q',x_LPT,fluid);
                s_LPT = CoolProp.PropsSI('S','T',T_LPT,'Q',x_LPT,fluid);
                ph_LPT = CoolProp.PropsSI('Phase','Q',x_LPT,'T',T_LPT,fluid);
            end
        elseif any(strcmp('T_LPT',varargin))  % If only T is available
            if ~isnumeric(varargin{find(strcmp('T_LPT',varargin))+1})
                error('T_LPT must be a number, not a %s.',class(varargin{find(strcmp('T_LPT',varargin))+1}))
            end
            T_LPT = varargin{find(strcmp('T_LPT',varargin))+1};
            x_LPT = x_LPT_default;
            P_LPT = CoolProp.PropsSI('P','T',T_LPT,'Q',x_LPT,fluid);
            rho_LPT = CoolProp.PropsSI('D','T',T_LPT,'Q',x_LPT,fluid);
            u_LPT = CoolProp.PropsSI('U','T',T_LPT,'Q',x_LPT,fluid);
            h_LPT = CoolProp.PropsSI('H','T',T_LPT,'Q',x_LPT,fluid);
            s_LPT = CoolProp.PropsSI('S','T',T_LPT,'Q',x_LPT,fluid);
            ph_LPT = CoolProp.PropsSI('Phase','P',P_LPT,'Q',x_LPT,fluid);
        else                                % If no variable for 4 is available
            T_LPT = T_LPT_default;
            x_LPT = x_LPT_default;
            P_LPT = CoolProp.PropsSI('P','T',T_LPT,'Q',x_LPT,fluid);
            rho_LPT = CoolProp.PropsSI('D','T',T_LPT,'Q',x_LPT,fluid);
            u_LPT = CoolProp.PropsSI('U','T',T_LPT,'Q',x_LPT,fluid);
            h_LPT = CoolProp.PropsSI('H','T',T_LPT,'Q',x_LPT,fluid);
            s_LPT = CoolProp.PropsSI('S','T',T_LPT,'Q',x_LPT,fluid);
            ph_LPT = CoolProp.PropsSI('Phase','P',P_LPT,'Q',x_LPT,fluid);
        end
    end  
    
%-------------------- Pump isoentropic efficiency ------------------------%    
    if any(strcmp('eta_p',varargin))
        if ~isnumeric(varargin{find(strcmp('eta_p',varargin))+1})
            error('eta_p must be a number, not a %s.',class(eta_p))
        end
        eta_p = varargin{find(strcmp('eta_p',varargin))+1};
    end

%--------------------- PRESSURE DROP IN THE HEATER -----------------------%
    if any(strcmp('DP_Hea',varargin))
        DeltaP_Hea = varargin{find(strcmp('DP_Hea',varargin))+1};
        if ~isnumeric(varargin{find(strcmp('DP_Hea',varargin))+1})
            error('DP_Hea must be a number, not a %s.',class(varargin{find(strcmp('DP_Hea',varargin))+1}))
        end
    else
        DeltaP_Hea = 0;
    end
    
%------------------------- Steady State analysis -----------------------%
% States:                                                               %
% 1 -> HPT outlet / evaporator inlet                                    %
% 2 -> evaporator outlet / turbine inlet                                %
% 3 -> turbine outlet / condenser inlet                                 %
% 4 -> condenser outlet / LPT inlet                                     %
% 5 -> LPT outlet / pump inlet                                          %
% 6 -> pump outlet / heater inlet                                       %
% 7 -> heater outlet / HTP inlet                                        %
%-----------------------------------------------------------------------%

% State at the outlet of the HPT
    P(1) = P_HPT;
    if ph_HPT == 0 % If the fluid is at the state of subcooled liquid in the HPT
        T(1) = T_HPT;
        h(1) = CoolProp.PropsSI('H','P',P(1),'T',T(1),fluid);
        s(1) = CoolProp.PropsSI('S','P',P(1),'T',T(1),fluid);
    elseif ph_HPT == 6 % If the fluid is at the state of liquid-vapour mixture in the HPT
        T(1) = CoolProp.PropsSI('T','P',P(1),'Q',0,fluid);
        h(1) = CoolProp.PropsSI('H','P',P(1),'Q',0,fluid);
        s(1) = CoolProp.PropsSI('S','P',P(1),'Q',0,fluid);
    else
        error('The fluid at the HPT must be either liquid or a liquid-vapour mixture.')
    end
    
% State at the evaporator outlet / turbine inlet
    P(2) = P(1) - DeltaP_Ev;
    T(2) = T_2;
    h(2) = CoolProp.PropsSI('H','P',P(2),'T',T(2),fluid);
    s(2) = CoolProp.PropsSI('S','P',P(2),'T',T(2),fluid);
    
% State at turbine outlet / condenser inlet
    P(3) = P_LPT;
    h_t = CoolProp.PropsSI('H','P',P(3),'S',s(2),fluid);
    h(3) = turbine('h_in',h(2),'h_s',h_t,'eta_t',eta_t);
    T(3) = CoolProp.PropsSI('T','P',P(3),'H',h(3),fluid);
    s(3) = CoolProp.PropsSI('S','P',P(3),'H',h(3),fluid);
    
% State at condenser outlet / pump inlet
    P(4) = P_LPT; %!!!
    x(4) = 0;
    T(4) = CoolProp.PropsSI('T','P',P(4),'Q',x(4),fluid);
    h(4) = CoolProp.PropsSI('H','P',P(4),'Q',x(4),fluid);
    s(4) = CoolProp.PropsSI('S','P',P(4),'Q',x(4),fluid);
    
% State at LPT outlet
    P(5) = P_LPT;
    x(5) = 0;
    T(5) = CoolProp.PropsSI('T','P',P(5),'Q',x(5),fluid);
    h(5) = CoolProp.PropsSI('H','P',P(5),'Q',x(5),fluid);
    s(5) = CoolProp.PropsSI('S','P',P(5),'Q',x(5),fluid);
    
% State at Pump outlet / Heater inlet
    P(6) = P_HPT; %!!!
    h_p = CoolProp.PropsSI('H','P',P(6),'S',s(5),fluid);
    h(6) = pump('h_in',h(5),'h_s',h_p,'eta_p',eta_p);
    T(6) = CoolProp.PropsSI('T','P',P(6),'H',h(6),fluid);
    s(6) = CoolProp.PropsSI('S','P',P(6),'H',h(6),fluid);
    
% State at Heater outlet / HPT inlet
    P(7) = P_HPT; %!!!
    x(7) = 0;
    T(7) = CoolProp.PropsSI('T','P',P(7),'Q',x(7),fluid);
    h(7) = CoolProp.PropsSI('H','P',P(7),'Q',x(7),fluid);
    s(7) = CoolProp.PropsSI('S','P',P(7),'Q',x(7),fluid);
    
% Energy transfer    
    q_Ev = h(2) - h(1);
%     q_C = h(3) - h(4)
    q_He = h(7) - h(6);
    w_t = h(2) - h(3);
    w_p = h(6) - h(5);

    switch Heat_source
        case 'HPs' % System with heat pump to provide heat at the evaporator/superheater
            
            % HEAT PUMP 1 - DISCHARGING PHASE
            % Hot source temperature (ORES system side of the heat exchanger)
            if T(2)>T_HPT % If there is superheating
                h_g_Ev = CoolProp.PropsSI('H','P',P(2),'Q',1,fluid);
                q_SH = h(2) - h_g_Ev;
                q_lg = q_Ev - q_SH;
                T_H1 = (q_lg/q_Ev)*T_HPT + (q_SH/q_Ev)*(T_HPT+T(2))/2;
            elseif T(2)==T_HPT % If saturated state
                T_H1 = T_HPT;
            else
                error('T_2 < T_sat')
            end

            COP_carnot_HP1 = 1/(1-T_amb/T_H1);    
            if COP_HP1_0 == 0
                COP_HP1 = COP_carnot_HP1;
            elseif COP_HP1_0 < COP_carnot_HP1
                COP_HP1 = COP_HP1_0;
            else
                error('COP_HP1 must be smaller than COP_carnot.')
            end
            w_c_HP1 = q_Ev/COP_HP1; % Work for reversible heat pump
            q_L_HP1 = q_Ev - w_c_HP1;      % Heat extracted from heat source

            % HEAT PUMP 2 - CHARGING PHASE
            T_H2 = (T(7)+T(6))/2;
%             T_H2 = T(7);
            COP_carnot_HP2 = 1/(1-T_amb/T_H2);
            if COP_HP2_0 == 0
                COP_HP2 = COP_carnot_HP2;
            elseif COP_HP2_0 < COP_carnot_HP2
                COP_HP2 = COP_HP2_0;
            else
                error('COP_HP2 must be smaller than COP_carnot.')
            end            
            w_c_HP2 = q_He/COP_HP2;
            q_L_HP2 = q_He - w_c_HP2;
            
            q_in = q_L_HP1 + q_L_HP2;
            w_net = w_t - w_p - w_c_HP1 - w_c_HP2;
            w_g = w_t - w_c_HP1;
%             w_g(1) = w_t(1);
            w_c = w_p + w_c_HP2;
        case 'Free'
            q_in = q_Ev + q_He;
            w_net = w_t - w_p;
            w_g = w_t;
            w_c = w_p;            
    end        
    
%---------------------------- Results display ----------------------------%
    if plotBool
        plotTS(fluid,'Iso_P',[P_HPT,P_LPT]);
        plotTS(fluid,'Process','Isobaric',P(1),s(1),s(2))
        plotTS(fluid,'Process','Isobaric',P(3),s(3),s(4))
        plotTS(fluid,'Process','Isobaric',P(7),s(6),s(7))
        plotTS(fluid,'Process','Expansion',P(2),h(2),P(3),eta_t)    
        plotTS(fluid,'Process','Pump_comp',P(5),h(5),P(6),eta_p)
        plot(s(1)./1000,T(1),'ob')
        plot(s(2)./1000,T(2),'ob')
        plot(s(3)./1000,T(3),'ob')
        plot(s(4)./1000,T(4),'ob')
        plot(s(5)./1000,T(5),'ob')
        plot(s(6)./1000,T(6),'ob')
        plot(s(7)./1000,T(7),'ob')     
%         plotPH(fluid)
%         hold all
%         plot(h(1)./1000,P(1)/1000,'o')
%         plot(h(2)./1000,P(2)/1000,'o')
%         plot(h(3)./1000,P(3)/1000,'o')
%         plot(h(4)./1000,P(4)/1000,'o')
%         plot(h(5)./1000,P(5)/1000,'o')
%         plot(h(6)./1000,P(6)/1000,'o')
%         plot(h(7)./1000,P(7)/1000,'o')       
    end
    
%=================== Full-cycle performance parameters ===================%
% Exergy
    ex = q_in*(1-T_amb/T_H)+w_p;
    %ex = (h(2) - h(5))- T_amb*(s(2)-s(5));
%     ex_HP = 
    eta_II = w_t/ex;
    
% Performance indexes
    eta_I = w_net/q_in;
%     eta_II = w_t/ex_HP;
    eta_RT = w_g/w_c;
    
    
% Cálculo alternativo da eficiência exergética
    % e_q_e = (h(2) - h(5))- T_amb*(s(2)-s(5));
    % eta_II = w_net/e_q_e;
    
end