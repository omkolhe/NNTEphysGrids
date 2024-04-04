function [] = rasterPlotPropColor(Spikes,prop1c,prop2sz,angleFlag)

if isempty(prop2sz)
    szFlag = 0;
else
    szFlag = 1;
end

y = length(Spikes(:,1));
x = length(Spikes(1,:));
M = load( 'myMap.mat' );
M.myMap = flip(M.myMap);

hold on;
for i=1:y
    x = find(Spikes(i,:)==1);
    y = i*ones(1,numel(x));
    if szFlag == 1
        scatter(x,y,prop2sz{1,i}, prop1c{1,i},"filled");
    else
        scatter(x,y,[], prop1c{1,i},"filled");
    end
end
colorbar;
if angleFlag==1 phasemap, else colormap("jet"), end
% phasemap;