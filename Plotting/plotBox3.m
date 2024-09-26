function plotBox3(data1,data2,data3)

% Making column vectors 
if size(data1,2) ~= 1
    data1 = data1';
end
if size(data2,2) ~= 1
    data2 = data2';
end
if size(data3,2) ~= 1
    data3 = data3';
end

n1 = size(data1,1);
n2 = size(data2,1);
n3 = size(data3,1);

maxN = max([n1,n2,n3]);

% Appending zeros
data1 = [data1;zeros(maxN-n1,1)];
data2 = [data2;zeros(maxN-n2,1)];
data3 = [data3;zeros(maxN-n3,1)];

customBoxplot([data1 data2 data3],'Scatter','on');

end