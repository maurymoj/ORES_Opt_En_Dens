function C_fluid = fluid_cost(fluid,mass,varargin)
% FLUID_COST Estimates fluid cost as a function of mass
%   fluid_cost()
% Inputs:
%   fluid
%   mass [kg]
% Data from Synquest - 2020 (http://synquestlabs.com/)

% Cost_index_ratio = (cost_index(2021)/cost_index(2020));  % Reference 2020
% because data was obtained in March 2020 2020 still has no average yearly
% CEPCI

if size(varargin)>0
    cost_source = varargin{1};
    switch cost_source
        case 'Synquest'
            
            mass = mass.*1000; % Converts mass to [g]

            switch fluid
                case 'R152a' % / 1,1-Difluoroethane
                    C_fluid = 1.232.*mass.^0.8671; % MATLAB Curve fitting
            %         C_fluid = 1.5519.*mass.^0.8489; % Excel trendline
                case 'R134a'
                    C_fluid = 4.892.*mass.^0.6351; % MATLAB Curve fitting
            %         C_fluid = 3.8255.*mass.^0.7043; % Excel trendline
                case 'R142b'
                    C_fluid = 12.3.*mass.^0.5754; % MATLAB Curve fitting
            %         C_fluid = 11.057.*mass.^0.6228; % Excel trendline
                case 'R365mfc'
                    C_fluid = 1.503.*mass.^0.7314; % MATLAB Curve fitting
            %         C_fluid = 0.9908.*mass.^0.8216; % Excel trendline
                case 'R236ea'
                    C_fluid = 1.806.*mass.^0.8802; % MATLAB Curve fitting
            %         C_fluid = 7.9118.*mass.^0.7004; % Excel trendline
                case 'R141b'
                    C_fluid = 18.64.*mass.^0.4952; % MATLAB Curve fitting
            %         C_fluid = 19.071.*mass.^0.5251; % Excel trendline
                otherwise
                    error('Fluid cost data not available.')     
            end
        case 'Alibaba'
            
            switch fluid
                case 'R152a' % / 1,1-Difluoroethane
                    C_fluid = mass.*2.96;
                case 'R134a'
                    C_fluid = mass.*4.78;
                case 'R142b'
                    C_fluid = mass.*3.06; 
                case 'R365mfc'
                    C_fluid = mass.*6.68;
%                 case 'R236ea'
%                     C_fluid = mass.;
                case 'R141b'
                    C_fluid = mass.*4.22;
                otherwise
                    error('Fluid cost data not available.')     
            end
       case 'Matheus'
            
            switch fluid
                case 'R152a' % / 1,1-Difluoroethane
                    C_fluid = mass.*7.01;
                case 'R134a'
                    C_fluid = mass.*4.49;
                case 'R142b'
                    C_fluid = mass.*4.41; 
                case 'R365mfc'
                    C_fluid = mass.*2;
%                 case 'R236ea'
%                     C_fluid = mass.;
                case 'R141b'
                    C_fluid = mass.*5;
                otherwise
                    error('Fluid cost data not available.')     
            end
        otherwise
            switch fluid
                case 'R152a' % / 1,1-Difluoroethane
                    C_fluid = mass.*2.96;
                case 'R134a'
                    C_fluid = mass.*4.78;
                case 'R142b'
                    C_fluid = mass.*3.06; 
                case 'R365mfc'
                    C_fluid = mass.*6.68;
%                 case 'R236ea'
%                     C_fluid = mass.;
                case 'R141b'
                    C_fluid = mass.*4.22;
                otherwise
                    error('Fluid cost data not available.')     
            end
    end
else
    switch fluid
        case 'R152a' % / 1,1-Difluoroethane
            C_fluid = mass.*2.96;
        case 'R134a'
            C_fluid = mass.*4.78;
        case 'R142b'
            C_fluid = mass.*3.06; 
        case 'R365mfc'
            C_fluid = mass.*6.68;
%                 case 'R236ea'
%                     C_fluid = mass.;
        case 'R141b'
            C_fluid = mass.*4.22;
        otherwise
            error('Fluid cost data not available.')     
    end
end