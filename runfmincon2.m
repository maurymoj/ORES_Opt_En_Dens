function [history,searchdir] = runfmincon2
 
% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];
searchdir = [];
 
% call optimization
x0 = [-1];
options = optimoptions(@fmincon,'OutputFcn',@outfun,... 
    'Display','iter','Algorithm','active-set');
xsol = fmincon(@objfun,x0,[],[],[],[],[],[],[],options);
 
 function stop = outfun(x,optimValues,state)
     stop = false;
 
     switch state
         case 'init'
             hold on
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x; x];
         % Concatenate current search direction with 
         % searchdir.
           searchdir = [searchdir;... 
                        optimValues.searchdirection'];
%            plot(x(1),'o');
%          % Label points with iteration number and add title.
%          % Add .15 to x(1) to separate label from plotted 'o'
%            text(x(1)+.15,... 
%                 num2str(optimValues.iteration));
%            title('Sequence of Points Computed by fmincon');
         case 'done'
             hold off
         otherwise
     end
 end
 
 function f = objfun(x)
     f= x(1).^4-4*x(1).^3+2.5*x(1).^2+2*x(1)-1;
 end
 
%  function [c, ceq] = confun(x)
%      % Nonlinear inequality constraints
%      c = [1.5 + x(1)*x(2) - x(1) - x(2);
%          -x(1)*x(2) - 10];
%      % Nonlinear equality constraints
%      ceq = [];
%  end
end