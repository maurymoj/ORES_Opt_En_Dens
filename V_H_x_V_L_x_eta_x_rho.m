clc;
fluids = {'R152a','R134a','R142b','R365mfc','R236ea','R141b'}; % ARTIGO
% Teste dos volumes dos tanques de alta e baixa pressão
K_V_H = 1.0:0.25:4;
K_V_L = 1.0:0.25:4;

max_RT = zeros(length(fluids),1);
max_rho = zeros(length(fluids),1);
min_RT = zeros(length(fluids),1);
min_rho = zeros(length(fluids),1);

for k = 5:5%length(fluids)
    fluid = fluids{k};
    eta_RT = zeros(length(K_V_H),length(K_V_L));
    rho_E = zeros(length(K_V_H),length(K_V_L));
    t_generation = zeros(length(K_V_H),length(K_V_L));
    erros = zeros(length(K_V_H),length(K_V_L));

    for i=1:length(K_V_H)
        for j = 1:length(K_V_L)
%             [K_V_H(i) K_V_L(j)]
            try
                [eta_RT(i,j),rho_E(i,j),t_generation(i,j)] = ORES_tr('K_V_HPT',K_V_H(i),'K_V_LPT',K_V_L(j),'fluid',fluid);%,'Plot',1)
            catch
                [i,j;K_V_H(i),K_V_L(j)] 
    %             disp('error')
%                 eta_RT(i,j) = NaN;
%                 rho_E(i,j) = NaN;
%                 t_generation(i,j) = 0;
                erros(i,j) = 1;
            end
        end    
    end

    eta_RT(t_generation<3600) = NaN;
    rho_E(t_generation<3600) = NaN;

    figure('color',[1 1 1])
    contour(K_V_L,K_V_H,eta_RT,'k','ShowText','on')
    hold on;
    contour(K_V_L,K_V_H,rho_E,'color',[0.7 0.7 0.7],'ShowText','on')
    grid on
    contour(K_V_L,K_V_H,t_generation==3600,[1,1],'r')
    xlabel('V_{LPT}')
    ylabel('V_{HPT}')
    legend('\eta_{RT}','\rho_E')
    % legend('\eta_{RT}','\rho_E','t_g == 3,600 s')
    title(fluid)
%     save(fluid)
    max_RT(k) = max(eta_RT(t_generation>=3600));
    max_rho(k) = max(rho_E(t_generation>=3600));
    min_RT(k) = min(eta_RT(t_generation>=3600));
    min_rho(k) = min(rho_E(t_generation>=3600));
end