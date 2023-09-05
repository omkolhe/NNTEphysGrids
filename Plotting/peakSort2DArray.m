function [xsorted,indexsort] = peakSort2DArray(x,direction,dimension)

xpeak = mean(x,dimension);

if dimension == 2
    [~,indexsort] = sort(xpeak,1,direction);
    xsorted = x(indexsort,:);
end
if dimension == 1
    [~,indexsort] = sort(xpeak,2,direction);
    xsorted = x(:,indexsort);
end

