function [chPointLinear] = chToLinear(chPointGrid,parameters)

% chPointGrid(1)  - columns
% chPointGrid(2) - rows

assert(( chPointGrid(1) <= parameters.cols && chPointGrid(2)<= parameters.rows), 'Check! Point out of actual grid coordinate' );

chPointLinear = (chPointGrid(2)-1)*parameters.cols + chPointGrid(1);
end

