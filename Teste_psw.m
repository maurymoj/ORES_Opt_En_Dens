fun = @dejong5fcn;
nvars = 2;
lb = [-50;-50];
ub = -lb;
rng default % For reproducibility
% options = optimoptions(@particleswarm,'OutputFcn',@pswplotranges);
options = optimoptions(@particleswarm,'OutputFcn',@pswplotmeanx);
[x,fval,exitflag] = particleswarm(fun,nvars,lb,ub,options)