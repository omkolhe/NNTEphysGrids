function [chPointGrid] = chToGrid(chPointLinear,parameters)

% chPointGrid(1)  - columns
% chPointGrid(2) - rows

assert(( chPointLinear <= parameters.rows*parameters.cols ), 'Check! Number of electrodes more than actual' );

chPointGrid(1) =  rem(chPointLinear,parameters.cols);
if chPointGrid(1) == 0, chPointGrid(1)=4; end
chPointGrid(2) = ceil(chPointLinear/parameters.cols);
end

