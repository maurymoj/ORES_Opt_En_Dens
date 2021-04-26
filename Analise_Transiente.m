rho_HPT_o = CoolProp.PropsSI('D','P',3500000,'Q',0.05,'R134a');
rho_HPT_f = CoolProp.PropsSI('D','P',3500000,'Q',0.95,'R134a');
Dm_min = 300;
V_HPT = Dm_min/(rho_HPT_o-rho_HPT_f);
rho_LPT_o = CoolProp.PropsSI('D','P',1000000,'Q',0.95,'R134a');
rho_LPT_f = CoolProp.PropsSI('D','P',1000000,'Q',0.05,'R134a');
V_LPT = Dm_min/(rho_LPT_f-rho_LPT_o);

%% ------------------------------------------------------------------------
% Teste do volume do tanque de alta pressão
K_V_H = 1:0.5:2.0;
% K_V_L = 0.5:0.2:1.5;

eta_RT = zeros(length(K_V_H),1);
rho_E = zeros(length(K_V_H),1);

for i=1:length(K_V_H)
%     [eta_RT(i),rho_E(i)] = ORES_tr('K_V_HPT',K_V_H(i),'Plot',1);
    [eta_RT(i),rho_E(i)] = ORES_tr('K_V_HPT',K_V_H(i));
%     for j = 1:length(K_V_L)
%        [eta_RT(i,j),rho_E(i,j)] = ORES_tr('K_V_HPT',K_V_H(i),'K_V_LPT',K_V_L(j));%,'Plot',1)
%     end    
end

%% Teste Dt
close all
% Dt = [5 10 25 50 100];
Dt = 5;
eta_RT = zeros(1,length(Dt));
rho_E = zeros(1,length(Dt));

for i=1:length(Dt)
    [eta_RT(i),rho_E(i),t_generation,~,CAPEX] = ORES_tr('K_V_HPT',1.5,'K_V_LPT',1.5,'fluid','R141b','Plot',1,'dt',Dt(i))
end
    
%% ------------------------------------------------------------------------
tic
clc;
fluids = {'R152a','R134a','R142b','R365mfc','R141b'}; % TESE
% fluids = {'R152a','R134a','R142b','R365mfc','R236ea','R141b'}; % TESE
% fluids = {'R152a','R134a','R142b','R365mfc','R141b'}; % ARTIGO
% Teste dos volumes dos tanques de alta e baixa pressão
K_V_H = 1.0:0.125:4;
K_V_L = 1.0:0.125:4;
% K_V_H = 1.0:1:4;
% K_V_L = 1.0:1:4;

max_RT = zeros(length(fluids),1);
max_rho = zeros(length(fluids),1);
min_RT = zeros(length(fluids),1);
min_rho = zeros(length(fluids),1);
% tic
for k = 1:length(fluids)
    fluid = fluids{k};
    eta_RT = zeros(length(K_V_H),length(K_V_L));
    rho_E = zeros(length(K_V_H),length(K_V_L));
    CAPEX = zeros(length(K_V_H),length(K_V_L));
    t_generation = zeros(length(K_V_H),length(K_V_L));
    erros = zeros(length(K_V_H),length(K_V_L));

    for i=1:length(K_V_H)
        for j = 1:length(K_V_L)
%             [K_V_H(i) K_V_L(j)]
            try
                [eta_RT(i,j),rho_E(i,j),t_generation(i,j),~,CAPEX(i,j)] = ORES_tr('K_V_HPT',K_V_H(i),'K_V_LPT',K_V_L(j),'fluid',fluid);%,'Plot',1)
            catch
%                 [i,j;K_V_H(i),K_V_L(j)] 
    %             disp('error')
                eta_RT(i,j) = NaN;
                rho_E(i,j) = NaN;
                CAPEX(i,j) = NaN;
%                 t_generation(i,j) = 0;
                erros(i,j) = 1;
            end
        end    
    end

    eta_RT(t_generation<3600) = NaN;
    rho_E(t_generation<3600) = NaN;
    CAPEX(t_generation<3600) = NaN;

    figure('color',[1 1 1])
    contour(K_V_L,K_V_H,eta_RT,'k','ShowText','on')
    hold on;
%     contour(K_V_L,K_V_H,rho_E,'color',[0.7 0.7 0.7],'ShowText','on')
    contour(K_V_L,K_V_H,CAPEX./1000,'color',[0.7 0.7 0.7],'ShowText','on')
    grid on
    contour(K_V_L,K_V_H,t_generation==3600,[1,1],'r')
    xlabel('K_{LPT}')
    xtickformat('%0.1f')
    ylabel('K_{HPT}')
    ytickformat('%0.1f')
    legend('\eta_{RT}','CAPEX')
    % legend('\eta_{RT}','\rho_E','t_g == 3,600 s')
    title(fluid)
    
    fig = gcf;
    fig.PaperOrientation = 'landscape';
    fig.PaperSize = [6 4.5];
%     name = strcat(fluid,'_',datestr(now,'yyyy_mm_dd'));
    name = strcat(fluid,'_','2020_03_31');
    save(name)
%     save(fluid)
    max_RT(k) = max(eta_RT(t_generation>=3600));
    max_rho(k) = max(rho_E(t_generation>=3600));
    min_RT(k) = min(eta_RT(t_generation>=3600));
    min_rho(k) = min(rho_E(t_generation>=3600));
end
% toc
%% Figura V x V x Custo x eta_RT

for k = 1:length(fluids)
    fluid = fluids{k};
%     name = strcat(fluid,'_',datestr(now,'yyyy_mm_dd'));
    name = strcat(fluid,'_2020_03_31.mat');
    load(name)

    figure('color',[1 1 1])
    contour(K_V_L,K_V_H,eta_RT,'k','ShowText','on')
    hold on;
    contour(K_V_L,K_V_H,CAPEX./1000,'color',[0.7 0.7 0.7],'ShowText','on')
    grid on
    xlabel('K_{LPT}')
    xtickformat('%0.1f')
    ylabel('K_{HPT}')
    ytickformat('%0.1f')
    legend('\eta_{RT}','CAPEX')
    title(fluid)
    
    fig = gcf;
    fig.PaperOrientation = 'landscape';
    fig.PaperSize = [6 4.5];
end

%% Figura V x V x rho_E x eta_RT

for k = 1:length(fluids)
    fluid = fluids{k};
%     name = strcat(fluid,'_',datestr(now,'yyyy_mm_dd'));
    name = strcat(fluid,'_2020_03_31.mat');
    load(name)

    figure('color',[1 1 1])
    contour(K_V_L,K_V_H,eta_RT,'k','ShowText','on')
    hold on;
    contour(K_V_L,K_V_H,rho_E,'color',[0.7 0.7 0.7],'ShowText','on')
%     contour(K_V_L,K_V_H,CAPEX./1000,'color',[0.7 0.7 0.7],'ShowText','on')
    grid on
%     contour(K_V_L,K_V_H,t_generation==3600,[1,1],'r')
    xlabel('K_{LPT}')
    xtickformat('%0.1f')
    ylabel('K_{HPT}')
    ytickformat('%0.1f')
    legend('\eta_{RT}','\rho_E')
    title(fluid)
    
    fig = gcf;
    fig.PaperOrientation = 'landscape';
    fig.PaperSize = [6 4.5];
end

%% Figura comparativa 
fluids = {'R152a','R134a','R142b','R365mfc','R141b'}; % TESE
fluids_std = {'R-152a','R-134a','R-142b','R-365mfc','R-141b'}; % TESE
% fluids = {'R152a','R134a','R142b','R365mfc','R236ea','R141b'}; % TESE
% fluids_std = {'R-152a','R-134a','R-142b','R-365mfc','R-236ea','R-141b'}; % TESE

% fluids = {'R152a','R134a','R142b','R365mfc','R141b'}; % ARTIGO
% fluids_std = {'R-152a','R-134a','R-142b','R-365mfc','R-141b'}; % ARTIGO
max_RT = zeros(length(fluids),1);
max_rho = zeros(length(fluids),1);
min_RT = zeros(length(fluids),1);
min_rho = zeros(length(fluids),1);
max_CAPEX = zeros(length(fluids),1);
min_CAPEX = zeros(length(fluids),1);

% figure('color',[1 1 1])
% applystyle2plot()
% hold on
marker = {'^','o','s','d','p','h'};

for k=1:length(fluids)
    k
    name = strcat(fluids{k},'_2020_03_31.mat');
    load(name);
    max_RT(k) = max(eta_RT(t_generation>=3600))
    max_rho(k) = max(rho_E(t_generation>=3600));
    max_CAPEX(k) = max(CAPEX(t_generation>=3600));
    min_RT(k) = min(eta_RT(t_generation>=3600));
    min_rho(k) = min(rho_E(t_generation>=3600));
    min_CAPEX(k) = min(CAPEX(t_generation>=3600));    
end
max_RT(1) = 0.4033;

% for k=1:length(fluids)
%     plot(max_rho(k),max_RT(k),marker{k},'MarkerEdgeColor','k','MarkerFaceColor','k')    
%     pause(1)
% end
% 
% for k=1:length(fluids)
%     plot(min_rho(k),min_RT(k),marker{k},'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])    
% end

% for k=1:length(fluids)
%     rectangle('Position',[min_rho(k) min_RT(k) (max_rho(k) - min_rho(k)) (max_RT(k) - min_RT(k))])
% end
% legend('','','','','','',fluids{1:length(fluids)},'location','NorthEast')
% grid on;
% ylim([0 1])
% xlabel('Energy density, \rho_E [Wh/L]')
% ylabel('Round-trip efficiency, \eta_{RT}')

% figure('color',[1 1 1])
% hold on
% grid on
% colors = {'r','g','b','y','m','c'};
% cor = 0.1:0.15:1;
% for k=1:length(fluids)
%     fill([min_rho(k) min_rho(k) max_rho(k) max_rho(k)],[min_RT(k) max_RT(k) max_RT(k) min_RT(k)],cor(k).*[1 1 1],'FaceAlpha',.5)
% end
% axis([0 2.5 0 1])
% ytickformat('%.1f')
% xtickformat('%.1f')
% xlabel('Energy density, \rho_E [Wh/L]')
% ylabel('Round-trip efficiency, \eta_{RT}')
% legend(fluids)


figure('color',[1 1 1])
hold on
grid on
applystyle2plot()

for k=1:length(fluids)
    plot(max_CAPEX(k)./1000,max_RT(k),marker{k},'MarkerEdgeColor','k','MarkerFaceColor','k')    
%     pause(1)
end

for k=1:length(fluids)
    plot(min_CAPEX(k)./1000,min_RT(k),marker{k},'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])    
end

% colors = {'r','g','b','y','m','c'};
cor = 0.1:0.15:1;
for k=1:length(fluids)
    fill([min_CAPEX(k) min_CAPEX(k) max_CAPEX(k) max_CAPEX(k)]./1000,[min_RT(k) max_RT(k) max_RT(k) min_RT(k)],cor(k).*[1 1 1],'FaceAlpha',.5)
end
legend(fluids{1:length(fluids)},'location','NorthEast')
% grid on;
ylim([0 1])
% axis([0 2.5 0 1])
ytickformat('%.1f')
xtickformat('%,1.0f')
xlabel('CAPEX [thousands of USD]')
ylabel('Round-trip efficiency, \eta_{RT}')
%% Comp V x V x eta x rho com fç fill
figure('color',[1 1 1])
hold on
grid on
max_RT = zeros(length(fluids),25);
for k=1:length(fluids)
    k
    name = strcat(fluids{k},'.mat');
    load(name);
%     if any(eta_RT(t_generation>=3600))
%         eta_RT(t_generation>=3600) = NaN;
%         rho_E(t_generation>=3600) = NaN;
        max_RT(k,:) = max(eta_RT);
        max_rho(k,:) = max(rho_E);
        min_RT(k,:) = min(eta_RT);
        min_rho(k,:) = min(rho_E);
%     end
end
% max_RT(1) = 0.4033;



%%
fluid = 'R141b';
K_V_H = 2;
K_V_L = 2;

P_L = CoolProp.PropsSI('P','T',298.15,'Q',0,fluid);
P_crit = CoolProp.Props1SI('Pcrit',fluid);

P_H_max = 0.8*P_crit;
P_H = 1.3*P_L:200000:P_H_max;

eta_RT = zeros(length(P_H),1);
rho_E = zeros(length(P_H),1);
t_generation = zeros(length(P_H),1);
erros = zeros(length(P_H),1);

for i=1:length(P_H)
%     try
        [eta_RT(i),rho_E(i),t_generation(i)] = ORES_tr('P_HPT',P_H(i),'K_V_HPT',K_V_H,'K_V_LPT',K_V_L,'fluid',fluid);%,'Plot',1)
%     catch
%         P_H(i)
%         disp('error')
%         eta_RT(i) = NaN;
%         rho_E(i) = NaN;
%         t_generation(i) = 0;
%         erros(i) = 1;
%     end
end

figure('color',[1 1 1])
plot(P_H./1000,eta_RT)
grid on;
hold on;
axis ([0 P_H_max./1000 0 1])
yyaxis right
plot(P_H./1000,rho_E)

%% Figura V x V x rho_E x eta_RT apresentação

for k = 1:length(fluids)
    fluid = fluids{k};
%     name = strcat(fluid,'_',datestr(now,'yyyy_mm_dd'));
    name = strcat(fluid,'_2020_03_31.mat');
    load(name)

    figure('color',[1 1 1])
    contour(K_V_L,K_V_H,eta_RT)
    xlabel('K_{LPT}')
    xtickformat('%0.1f')
    ylabel('K_{HPT}')
    ytickformat('%0.1f')
    title(fluid)
    
    fig = gcf;
    fig.PaperOrientation = 'landscape';
    fig.PaperSize = [6 4.5];
    
    figure('color',[1 1 1])
    contour(K_V_L,K_V_H,rho_E)
%     contour(K_V_L,K_V_H,CAPEX./1000,'color',[0.7 0.7 0.7],'ShowText','on')
%     grid on
%     contour(K_V_L,K_V_H,t_generation==3600,[1,1],'r')
    xlabel('K_{LPT}')
    xtickformat('%0.1f')
    ylabel('K_{HPT}')
    ytickformat('%0.1f')
%     legend('\eta_{RT}','\rho_E')
    title(fluid)
    
    fig = gcf;
    fig.PaperOrientation = 'landscape';
    fig.PaperSize = [6 4.5];
end