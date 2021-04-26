function [varargout] = plotTS( fluid, varargin )
%plotTS Plots the T-s diagram for a fluid. Hold plot already on.
%   plotTS(Fluid,varargin)
%   if there is already a defined lower limit for the temperature it is
%   given as the second input, otherwise the lower limit is set to the 
%   temperature at the triple point
% DESCREVER VARIAÇÕES
if any(strcmp('Process',varargin))
    proc = varargin{find(strcmp('Process',varargin))+1};    
    if ~ischar(proc)
        error('Process must be a string, not a %s.',class(proc))
    end
    
    n_subd = 100;
    switch proc
        case 'Isobaric'
            Iso_P = varargin{find(strcmp('Process',varargin))+2};
            s1 = varargin{find(strcmp('Process',varargin))+3};
            s2 = varargin{find(strcmp('Process',varargin))+4};
            if ~isnumeric(Iso_P)
                error('Iso_P must be a number, not a %s.',class(Iso_P))
            elseif ~isnumeric(s1)
                error('s_1 must be a number, not a %s.',class(s1))
            elseif ~isnumeric(s2)
                error('s_2 must be a number, not a %s.',class(s2))
            end

            s_Iso_P = s1:(s2-s1)/(n_subd-1):s2;
            T_Iso_P = zeros(1,n_subd);
            for j=1:n_subd
                T_Iso_P(j) = CoolProp.PropsSI('T','P',Iso_P,'S',s_Iso_P(j),fluid);
            end
            
            if any(strcmp('plotColor',varargin))
                plotColor = varargin{find(strcmp('plotColor',varargin))+1};
                if ~ischar(plotColor) && ~isvector(plotColor)
                    error('plotColor must be a string or a vector, not a %s.',class(plotColor))
                elseif isnumeric(plotColor) && isvector(plotColor) && length(plotColor)~=3
                    error('plotColor must have 3 elements.')
                elseif ischar(plotColor)
                    handle=plot(s_Iso_P./1000,T_Iso_P,plotColor);
                    if nargout == 1
                        varargout{1} = handle;
                    end
                elseif isvector(plotColor)
                    handle=plot(s_Iso_P./1000,T_Iso_P,'color',plotColor);
                    if nargout == 1
                        varargout{1} = handle;
                    end
                end
            else
                plotColor = 'b';
                handle=plot(s_Iso_P./1000,T_Iso_P,plotColor);
                if nargout == 1
                    varargout{1} = handle;
                end
            end            
        case 'Expansion'
            P1 = varargin{find(strcmp('Process',varargin))+2};
            h1 = varargin{find(strcmp('Process',varargin))+3};
            P2 = varargin{find(strcmp('Process',varargin))+4};            
            eta_t = varargin{find(strcmp('Process',varargin))+5};
            if ~isnumeric(P1)
                error('P_1 must be a number, not a %s.',class(P1))
            elseif ~isnumeric(h1)
                error('h_1 must be a number, not a %s.',class(h1))
            elseif ~isnumeric(P2)
                error('P_2 must be a number, not a %s.',class(P2))
            elseif ~isnumeric(eta_t)
                error('eta must be a number, not a %s.',class(eta_t))
            end

            s_exp = zeros(1,n_subd);
            T_exp = zeros(1,n_subd);
            h_exp = zeros(1,n_subd);
            P_exp = P1:-(P1-P2)/(n_subd-1):P2;
            h_exp(1) = h1;          
            T_exp(1) = CoolProp.PropsSI('T','P',P_exp(1),'H',h_exp(1),fluid);
            s_exp(1) = CoolProp.PropsSI('S','P',P_exp(1),'H',h_exp(1),fluid);            
            
            for j=2:n_subd
                h_t = CoolProp.PropsSI('H','P',P_exp(j),'S',s_exp(j-1),fluid);
                h_exp(j) = turbine('h_in',h_exp(j-1),'h_s',h_t,'eta_t',eta_t);
                T_exp(j) = CoolProp.PropsSI('T','P',P_exp(j),'H',h_exp(j),fluid);
                s_exp(j) = CoolProp.PropsSI('S','P',P_exp(j),'H',h_exp(j),fluid);
            end
            
            if any(strcmp('plotColor',varargin))
                plotColor = varargin{find(strcmp('plotColor',varargin))+1};
                if ~ischar(plotColor) && ~isvector(plotColor)
                    error('plotColor must be a string or a vector, not a %s.',class(plotColor))
                elseif isnumeric(plotColor) && isvector(plotColor) && length(plotColor)~=3
                    error('plotColor must have 3 elements.')
                elseif ischar(plotColor)
                    handle=plot(s_exp./1000,T_exp,plotColor);
                    if nargout == 1
                        varargout{1} = handle;
                    end
                elseif isvector(plotColor)
                    handle=plot(s_exp./1000,T_exp,'color',plotColor);
                    if nargout == 1
                        varargout{1} = handle;
                    end
                end
            else
                plotColor = 'b';
                handle=plot(s_exp./1000,T_exp,plotColor);
                if nargout == 1
                    varargout{1} = handle;
                end
            end
        case 'Pump_comp'
% CONFERIR !!!!!!!!!!!!!!!!!!!!!!!!!!
            P1 = varargin{find(strcmp('Process',varargin))+2};
            h1 = varargin{find(strcmp('Process',varargin))+3};
            P2 = varargin{find(strcmp('Process',varargin))+4};            
            eta_p = varargin{find(strcmp('Process',varargin))+5};
            if ~isnumeric(P1)
                error('P_1 must be a number, not a %s.',class(P1))
            elseif ~isnumeric(h1)
                error('h_1 must be a number, not a %s.',class(h1))
            elseif ~isnumeric(P2)
                error('P_2 must be a number, not a %s.',class(P2))
            elseif ~isnumeric(eta_p)
                error('eta must be a number, not a %s.',class(eta_p))
            end

            s_pump = zeros(1,n_subd);
            T_pump = zeros(1,n_subd);
            h_pump = zeros(1,n_subd);
            P_pump = P1:(P2-P1)/(n_subd-1):P2;
            h_pump(1) = h1;
            % if phase is saturated liquid A elseif subcooled B
            s_pump(1) = CoolProp.PropsSI('S','P',P_pump(1),'H',h_pump(1),fluid);
            T_pump(1) = CoolProp.PropsSI('T','P',P_pump(1),'H',h_pump(1),fluid);
            
            for j=2:n_subd
                h_s = CoolProp.PropsSI('H','P',P_pump(j),'S',s_pump(j-1),fluid);
                h_pump(j) = pump('h_in',h_pump(j-1),'h_s',h_s,'eta_p',eta_p);
                T_pump(j) = CoolProp.PropsSI('T','P',P_pump(j),'H',h_pump(j),fluid);
                s_pump(j) = CoolProp.PropsSI('S','P',P_pump(j),'H',h_pump(j),fluid);
            end
            
            if any(strcmp('plotColor',varargin))
                plotColor = varargin{find(strcmp('plotColor',varargin))+1};
                if ~ischar(plotColor) && ~isvector(plotColor)
                    error('plotColor must be a string or a vector, not a %s.',class(plotColor))
                elseif isnumeric(plotColor) && isvector(plotColor) && length(plotColor)~=3
                    error('plotColor must have 3 elements.')
                elseif ischar(plotColor)
                    handle=plot(s_pump./1000,T_pump,plotColor);
                    if nargout == 1
                        varargout{1} = handle;
                    end
                elseif isvector(plotColor)
                    handle=plot(s_pump./1000,T_pump,'color',plotColor);
                    if nargout == 1
                        varargout{1} = handle;
                    end
                end
            else
                plotColor = 'b';
                handle=plot(s_pump./1000,T_pump,plotColor);
                if nargout == 1
                    varargout{1} = handle;
                end
            end
    end
%     n_subd = 100;
%     switch proc
%         case 'Isobaric'
%             Iso_P = varargin{find(strcmp('Process',varargin))+2};
%             s1 = varargin{find(strcmp('Process',varargin))+3};
%             s2 = varargin{find(strcmp('Process',varargin))+4};
%             if ~isnumeric(Iso_P)
%                 error('Iso_P must be a number, not a %s.',class(Iso_P))
%             elseif ~isnumeric(s1)
%                 error('s_1 must be a number, not a %s.',class(s1))
%             elseif ~isnumeric(s2)
%                 error('s_2 must be a number, not a %s.',class(s2))
%             end
%             
%             s_Iso_P = s1:(s2-s1)/(n_subd-1):s2;
%             T_Iso_P = zeros(1,n_subd);
%             for j=1:n_subd
%                 T_Iso_P(j) = CoolProp.PropsSI('T','P',Iso_P,'S',s_Iso_P(j),fluid);
%             end
%             plot(s_Iso_P./1000,T_Iso_P,'b')
%         case 'Expansion'
%             P1 = varargin{find(strcmp('Process',varargin))+2};
%             h1 = varargin{find(strcmp('Process',varargin))+3};
%             P2 = varargin{find(strcmp('Process',varargin))+4};            
%             eta_t = varargin{find(strcmp('Process',varargin))+5};
%             if ~isnumeric(P1)
%                 error('P_1 must be a number, not a %s.',class(P1))
%             elseif ~isnumeric(h1)
%                 error('h_1 must be a number, not a %s.',class(h1))
%             elseif ~isnumeric(P2)
%                 error('P_2 must be a number, not a %s.',class(P2))
%             elseif ~isnumeric(eta_t)
%                 error('eta must be a number, not a %s.',class(eta_t))
%             end
%             
%             s_exp = zeros(1,n_subd);
%             T_exp = zeros(1,n_subd);
%             h_exp = zeros(1,n_subd);
%             P_exp = P1:-(P1-P2)/(n_subd-1):P2;
%             h_exp(1) = h1;          
%             T_exp(1) = CoolProp.PropsSI('T','P',P_exp(1),'H',h_exp(1),fluid);
%             s_exp(1) = CoolProp.PropsSI('S','P',P_exp(1),'H',h_exp(1),fluid);            
%             
%             for j=2:n_subd
%                 h_t = CoolProp.PropsSI('H','P',P_exp(j),'S',s_exp(j-1),fluid);
%                 h_exp(j) = turbine('h_in',h_exp(j-1),'h_s',h_t,'eta_t',eta_t);
%                 T_exp(j) = CoolProp.PropsSI('T','P',P_exp(j),'H',h_exp(j),fluid);
%                 s_exp(j) = CoolProp.PropsSI('S','P',P_exp(j),'H',h_exp(j),fluid);
%             end
%             plot(s_exp./1000,T_exp,'b')
%         case 'Pump_comp'
%             P1 = varargin{find(strcmp('Process',varargin))+2};
%             h1 = varargin{find(strcmp('Process',varargin))+3};
%             P2 = varargin{find(strcmp('Process',varargin))+4};            
%             eta_p = varargin{find(strcmp('Process',varargin))+5};
%             if ~isnumeric(P1)
%                 error('P_1 must be a number, not a %s.',class(P1))
%             elseif ~isnumeric(h1)
%                 error('h_1 must be a number, not a %s.',class(h1))
%             elseif ~isnumeric(P2)
%                 error('P_2 must be a number, not a %s.',class(P2))
%             elseif ~isnumeric(eta_p)
%                 error('eta must be a number, not a %s.',class(eta_p))
%             end
%             
%             s_exp = zeros(1,n_subd);
%             T_exp = zeros(1,n_subd);
%             h_exp = zeros(1,n_subd);
%             P_exp = P1:(P2-P1)/(n_subd-1):P2;
%             h_exp(1) = h1;
%             if phase is saturated liquid A elseif subcooled B
%             s_exp(1) = CoolProp.PropsSI('S','P',P_exp(1),'H',h_exp(1),fluid);
%             T_exp(1) = CoolProp.PropsSI('T','P',P_exp(1),'H',h_exp(1),fluid);
%             
%             for j=2:n_subd
%                 h_p = CoolProp.PropsSI('H','P',P_exp(j),'S',s_exp(j-1),fluid);
%                 h_exp(j) = pump('h_in',h_exp(j-1),'h_s',h_p,'eta_p',eta_p);
%                 T_exp(j) = CoolProp.PropsSI('T','P',P_exp(j),'H',h_exp(j),fluid);
%                 s_exp(j) = CoolProp.PropsSI('S','P',P_exp(j),'H',h_exp(j),fluid);
%             end
%             plot(s_exp./1000,T_exp,'r')
%     end
else
    T_crit = CoolProp.Props1SI('Tcrit',fluid);

    if any(strcmp('T_min',varargin))
        T_min = varargin{find(strcmp('T_min',varargin))+1};
        if ~isnumeric(P_min)
            error('T_min must be a number, not a %s.',class(T_min))
        end
    else
        T_min = CoolProp.Props1SI(fluid,'T_TRIPLE');
    end
    
    T = T_min:T_crit;
    s = zeros(2*length(T),1);
    for i=1:length(T)
        s(i)=CoolProp.PropsSI('S','T',T(i),'Q',0,fluid);
        s(i+length(T))=CoolProp.PropsSI('S','T',T(length(T)+1-i),'Q',1,fluid);
    end

    n=length(T);
    for i=1:n
        T(n+i)=T(n+1-i);
    end

    handle=figure('color',[1 1 1]);
    if nargout == 1
        varargout{1} = handle;
    end
    plot(s/1000,T,'k')
    xlabel('Entropy [kJ/{kg-K}]')
    ylabel('Temperature [K]')
    grid on;
    hold on;
end
%----------------------- Inclusao de iso-linhas --------------------------%
n_subd = 100;
if any(strcmp('Iso_P',varargin))
    Iso_P = varargin{find(strcmp('Iso_P',varargin))+1};
    if ~isnumeric(Iso_P)
        error('Iso_P must be a number, not a %s.',class(Iso_P))
    elseif ~isvector(Iso_P) 
        error('Iso_P must be a vector.')
    end
    s_Iso_P = min(s):(max(s)-min(s))/(n_subd-1):max(s);
    T_Iso_P = zeros(length(Iso_P),n_subd);
    for i=1:length(Iso_P)
        %i
        for j=1:n_subd
            %j
            T_Iso_P(i,j) = CoolProp.PropsSI('T','P',Iso_P(i),'S',s_Iso_P(j),fluid);
        end
        handle=plot(s_Iso_P./1000,T_Iso_P(i,:),'--','color',[0.6 0.6 0.6]);
        if nargout == 1
            varargout{1} = handle;
        end
    end
end
if any(strcmp('Iso_v',varargin))
    Iso_v = varargin{find(strcmp('Iso_v',varargin))+1};
    if ~isnumeric(Iso_v)
        error('Iso_v must be a number, not a %s.',class(Iso_v))
    elseif ~isvector(Iso_v) 
        error('Iso_v must be a vector.')
    end
    s_Iso_v = min(s):(max(s)-min(s))/(n_subd-1):max(s);
    T_Iso_v = zeros(length(Iso_v),n_subd);
    for i=1:length(Iso_v)
        %i
        for j=1:n_subd
            %j
            T_Iso_v(i,j) = CoolProp.PropsSI('T','S',s_Iso_v(j),'D',1/Iso_v(i),fluid);
        end
        handle=plot(s_Iso_h./1000,T_Iso_v(i,:),'--','color',[0.6 0.6 0.6]);
        if nargout == 1
            varargout{1} = handle;
        end
    end
end

%-------------------------------------------------------------------------%
end