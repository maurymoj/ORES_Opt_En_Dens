% otim_ORES_st_DTSHxP
%   Maximização da eta_RT tendo como variável de decisão a pressão P_HPT
P = [1500:250:3500]'; % P´[kPa]
X = zeros(length(P),2);
eta_RT = zeros(length(P),1);

for i=1:length(P)
    X0 = [P(i),20];
    LB = [P(i),0.01];
    UB = [P(i),40];

    [X(i,:),eta_RT(i)] = fmincon(@otim_ORES_st,X0,[],[],[],[],LB,UB);    
end

figure('color',[1 1 1])
grid on;
hold all;
%plot(X(:,1),eta_RT)
[hAx,hLine1,hLine2] = plotyy(X(:,1),-eta_RT,X(:,1),X(:,2));
xlabel('Pressure [kPa]')

ylabel(hAx(1),'\eta_{RT}') % left y-axis 
ylabel(hAx(2),'\Delta T_{SH}^*') % right y-axis

%% "Case study"
ind = 9;
ORES_st('R134a','Plot',1,'P_HPT',X(ind,1)*1000,'DT_SH',X(ind,2))