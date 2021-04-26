function [handle] = plotPH( fluid, varargin )
%PLOTPH Plots the P-h diagram for 'Fluid', 'hold' already set to 'on'.
%   if there is already a defined lower limit for the pressure it is
%   given as the second input, otherwise the lower limit is set to the 
%   pressure at the triple point
P_crit = CoolProp.Props1SI('Pcrit',fluid);

if any(strcmp('P_min',varargin))
    P_min = varargin{find(strcmp('P_min',varargin))+1};
    if ~isnumeric(P_min)
        error('P_min must be a number, not a %s.',class(P_min))
    end
else
    P_min = CoolProp.Props1SI(fluid,'P_TRIPLE');
end

P = P_min:1000:P_crit-30000;
h = zeros(2*length(P),1);
for i=1:length(P)
    h(i)=CoolProp.PropsSI('H','P',P(i),'Q',0,fluid);
    h(i+length(P))=CoolProp.PropsSI('H','P',P(length(P)+1-i),'Q',1,fluid);
end

n=length(P);
for i=1:n
    P(n+i)=P(n+1-i);
end

handle=figure('color',[1 1 1]);
plot(h./1000,P./1000,'k')
xlabel('Enthalpy (kJ/kg)')
ylabel('Pressure [kPa]')
set(gca, 'YScale', 'log')
grid on;
hold on;

% %----------------------- Inclusão de iso-linhas --------------------------%
% n_subd = 100;
% if any(strcmp('Iso_T',varargin))
%     Iso_T = varargin{find(strcmp('Iso_T',varargin))+1};
%     if ~isnumeric(Iso_T)
%         error('Iso_T must be a number, not a %s.',class(Iso_T))
%     elseif ~isvector(Iso_T) 
%         error('Iso_T must be a vector.')
%     end
%     h_Iso_T = min(h):(max(h)-min(h))/(n_subd-1):max(h);
%     P_Iso_T = zeros(length(Iso_T),n_subd);
%     for i=1:length(Iso_T)
%         i
%         for j=1:n_subd
%             ERRO; PAR DE ENTRADA h-T AINDA NÃO ESTÁ HABILITADO
%             P_Iso_T(i,j) = CoolProp.PropsSI('P','T',Iso_T(i),'H',h_Iso_T(j),fluid);
%         end
%         plot(h_Iso_T(i,:)./1000,P_Iso_T(i,:)./1000,'--')
%     end
% end

end