PreCueBaselinePGD = mean(PGD.PGDHitsBaseline(:,1000:1500),2);
PreCueBaseline = mean(PreCueBaselinePGD);
PerPreCueBaselinePGD = 100*(mean(PGD.PGDHitsBaseline(:,1000:1500),2)-PreCueBaseline)/PreCueBaseline;
PerPostCueBaselinePGD = 100*((mean(PGD.PGDHitsBaseline(:,1500:2000),2) - PreCueBaselinePGD)./PreCueBaselinePGD);
% figure,plotBox2(PerPreCueBaselinePGD,PerPostCueBaselinePGD);

PreCueCoolingPGD = mean(PGD.PGDHitsCooling(:,1000:1500),2);
PreCueBaseline = mean(PreCueCoolingPGD);
PerPreCueCoolingPGD = 100*(mean(PGD.PGDHitsCooling(:,1000:1500),2)-PreCueBaseline)/PreCueBaseline;
PerPostCueCoolingPGD = 100*((mean(PGD.PGDHitsCooling(:,1500:2000),2) - PreCueCoolingPGD)./PreCueCoolingPGD);
% figure,plotBox2(PerPreCueCoolingPGD,PerPostCueCoolingPGD);

% Comparing changes 
ranksum(PerPreCueBaselinePGD,PerPostCueBaselinePGD)
ranksum(PerPreCueCoolingPGD,PerPostCueCoolingPGD)

% Comparing both pres and posts
ranksum(PerPreCueBaselinePGD,PerPreCueCoolingPGD)
ranksum(PerPostCueBaselinePGD,PerPostCueCoolingPGD)


figure();
connectedBoxPlot(PerPreCueBaselinePGD,PerPostCueBaselinePGD);
ylim([-35 70]);box off;
figure();
connectedBoxPlot(PerPreCueCoolingPGD,PerPostCueCoolingPGD);
ylim([-35 70]);box off;
