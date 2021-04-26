fluids = {'R152a','R134a','R142b','R365mfc','R236ea','R141b'}; % ARTIGO

for k=1:length(fluids)
    plotTS(fluids{k})
    title(fluids{k})
end

% Conclus�es: 
% R152a e R134a - molhados
% R141b e R142b - ~isoentr�picos
% R365mfc e R236ea - secos