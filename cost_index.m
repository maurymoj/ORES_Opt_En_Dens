function cost_index = cost_index(year,varargin)
% COST_INDEX returns the cost index based on year
%   cost_index(YEAR) = returns the cost index for YEAR on CEPCI base
%   cost_index(YEAR,BASE) = return the cost index from the selected base
%   (currently only two bases are available - Chemical Engineering Plant 
%   Cost Index and Marshall and Swift Process Industry index)
if ~isempty(varargin)
    base = varargin{1};
else
    base = '';
end

% Data obtained from Turton et al., 2018 and https://www.chemengonline.com/2019-cepci-updates-january-prelim-and-december-2018-final/
CEPCI_index_array = [382 387 390 391 394 394 396 402 444 468 500 525 575 521 551 586 585 567 576 557 542 568 603.1 607.5];
Marshall_index_array = [1036 1053 1062 1062 1070 1095 1096 1113 1133 1218 1275 1354 1393 1487 1447 1477 1537 1553 1567 1598 1582];
switch base
    case 'CEPCI'
        if year >= 1996 & year <= 2019
            cost_index = CEPCI_index_array(year - 1995);
        else
            error('Index not available for the year requested.')
        end
    case 'Marshall and Swift'
        if year >= 1996 & year <= 2016
            cost_index = Marshall_index_array(year - 1995);
        else
            error('Index not available for the year requested.')
        end
    otherwise
        if year >= 1996 & year <= 2019
            cost_index = CEPCI_index_array(year - 1995);
        else
            error('Index not available for the year requested.')
        end
end

end