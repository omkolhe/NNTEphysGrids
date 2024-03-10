function plotBox2(data1,data2)

% Making column vectors 
if size(data1,2) ~= 1
    data1 = data1';
end
if size(data2,2) ~= 1
    data2 = data2';
end

n1 = size(data1,1);
n2 = size(data2,1);
% Appending zeros
if n1>n2
    data2 = [data2;zeros(n1-n2,1)];
else
    data1 = [data1;zeros(n2-n1,1)];
end 

customBoxplot([data1 data2],'Scatter','on');

end