function [history] = runpswarm
% , searchdir
% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];
% searchdir = [];

fun = @dejong5fcn;
nvars = 2;
lb = [-50;-50];
ub = -lb;

% call optimization
rng default % For reproducibility
options = optimoptions(@particleswarm,'OutputFcn',@outfun,... 
    'Display','iter');
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
