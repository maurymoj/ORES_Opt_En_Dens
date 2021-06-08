function [history] = runga_ORES

% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];
% searchdir = [];

fluid = 'R365mfc';
P_min = CoolProp.PropsSI('P','T',298.15,'Q',0,fluid);
P_max = 0.95*CoolProp.Props1SI(fluid,'Pcrit');

fun = @ORES_tr_opt;
nvars = 3;
lb = [1.2;1.2;1.5*P_min];
ub = [4;4;P_max];


% call optimization
rng default % For reproducibility
options = optimoptions(@ga,'OutputFcn',@outfun,... 
    'Display','iter');
[x,fval] = ga(fun,nvars,[],[],[],[],lb,ub,@mycon,options)

function [state,options,optchanged] = outfun(options,state,flag)
     
optchanged = false; 
     switch flag
         case 'init'
             hold on
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           history.fval = [history.fval, state.Score(:)];
           history.x = [history.x, state.Population(:)];
         % Concatenate current search direction with 
         % searchdir.
%            searchdir = [searchdir,... 
%                         optimValues.searchdirection'];
%            plot(x(1),x(2),'o');
         % Label points with iteration number and add title.
         % Add .15 to x(1) to separate label from plotted 'o'
%            text(x(1)+.15,x(2),... 
%                 num2str(optimValues.iteration));
%            title('Sequence of Points Computed by fmincon');
            if state.Generation == 25
                options.CrossoverFraction = 0.8;
                optchanged = true;
            end
         case 'done'
             hold off
         otherwise
     end
end


end