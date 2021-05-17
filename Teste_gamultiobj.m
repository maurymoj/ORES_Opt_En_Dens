%%
% fitfun = @(x)[-x^2,x^2-2];
fitfun = @(x)[sin(x),cos(x)];
nvars = 1;

lb = -2;
ub = 2;

options = optimoptions('gamultiobj','PlotFcn',@gaplotpareto);
% x = gamultiobj(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options)
[x,fval,exitflag,output,population,scores] = gamultiobj(fitfun,nvars,[],[],[],[],lb,ub,@mycon,options);

plot(sin(x),cos(x),'ko')
xlabel('f1')
ylabel('f2')
title('Pareto Front')


%%


rng default

% fitnessfcn = @(x)[norm(x)^2,0.5*norm(x(:)-[2;-1])^2+2];
nvars = 2;

% A = [1,1];
% b = 1/2;

% lb = 0;
% ub = 2*pi;

[x,fval,exitflag,output,population,scores] = gamultiobj(fitnessfcn,nvars,...
    [],[],[],[],lb,ub);

plot(x(:,1),x(:,2),'ko')
% t = linspace(-1/2,2);
% y = 1/2 - t;
% hold on
% plot(t,y,'b--')
xlabel('x(1)')
ylabel('x(2)')
title('Pareto Points in Parameter Space')

%% 
rng default
fitnessfcn = @(x)[sin(x),cos(x)];
nvars = 1;
lb = 0;
up = 2*pi;
x = gamultiobj(fitnessfcn,nvars,[],[],[],[],lb,ub)

plot(sin(x),cos(x),'r*')
xlabel('sin(x)')
ylabel('cos(x)')
title('Pareto Front')
legend('Pareto front')