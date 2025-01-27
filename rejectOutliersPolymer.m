[~,TFBaselineReject] = rmoutliers(ISIBaseline);
[~,TFPolymerReject] = rmoutliers(ISIPolymer);

ISIBaselineReject = ISIBaseline(~(TFPolymerReject | TFBaselineReject));
ISIPolymerReject = ISIPolymer(~(TFPolymerReject | TFBaselineReject));

temp = [ISIBaselineReject';ISIPolymerReject'];
templabels = cellstr([repmat('Baseline',size(ISIBaselineReject,2),1);repmat('Polymer ',size(ISIPolymerReject,2),1)]);
colors = [166/255 14/255 90/255;53/255 189/255 206/255];
figure,violinplot(temp,templabels,'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
ylabel('ISI (s)'); title('Polymer nPBDF');
% xtix = {'Baseline','Cooled','Recovery'}; xtixloc = [1 2 3]; set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);set(gca,'TickDir','out','fontsize',14');
set(gca,'TickDir','out','fontsize',14');box off;

ranksum(a,b)
signrank(ISIBaselineReject,ISIPolymerReject)

figure();
connectedBoxPlot(ISIBaselineReject,ISIPolymerReject);
ylim([-35 70]);box off;

