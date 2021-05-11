
rng default

fitnessfcn = @(x)[norm(x)^2,0.5*norm(x(:)-[2;-1])^2+2];
nvars = 2;

A = [1,1];
b = 1/2;

lb = 0;
ub = 2*pi;


[x,fval,exitflag,output,population,scores] = gamultiobj(fitnessfcn,nvars,...
    [],[],[],[],lb,ub);

plot(x(:,1),x(:,2),'ko')
xlabel('x(1)')
ylabel('x(2)')
title('Pareto Points in Parameter Space')