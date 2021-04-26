fluids = {'R152a','R134a','R142b','R365mfc','R236ea','R141b'}; % ARTIGO

for k=1:length(fluids)
    plotTS(fluids{k})
    title(fluids{k})
end

% Conclusões: 
% R152a e R134a - molhados
% R141b e R142b - ~isoentrópicos
% R365mfc e R236ea - secos