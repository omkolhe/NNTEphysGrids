function connectedBoxPlot(data1,data2)

% Making column vectors 
if size(data1,2) ~= 1
    data1 = data1';
end
if size(data2,2) ~= 1
    data2 = data2';
end

assert((size(data1,1)==size(data2,1)),'Variable size does not match')

n = size(data1,1);

boxplot([data1,data2]);
hold on;

x1 = ones(n,1);
x2 = 2*ones(n,1);
line([x1 x2]',[data1 data2]','Color',[0.7 0.7 0.7],'LineWidth',1);

end