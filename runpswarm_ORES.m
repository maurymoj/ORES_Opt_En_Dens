function [history] = runpswarm_ORES

% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];
% searchdir = [];

fluid = 'R141b';
P_min = CoolProp.PropsSI('P','T',298.15,'Q',0,fluid);
P_max = 0.95*CoolProp.Props1SI(fluid,'Pcrit');

fun = @ORES_tr_opt;
nvars = 3;
lb = [1;1;P_min];
ub = [4;4;P_max];

% call optimization
rng default % For reproducibility
options = optimoptions(@particleswarm,'OutputFcn',@outfun,... 
    'Display','iter','MaxIterations',50);
[x,fval,exitflag] = particleswarm(fun,nvars,lb,ub,options)

function stop = outfun(optimValues,state)
     stop = false;
 
     switch state
         case 'init'
             hold on
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           history.fval = [history.fval, optimValues.swarmfvals];
           history.x = [history.x, optimValues.swarm(:)];
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
         case 'done'
             hold off
         otherwise
     end
end


end