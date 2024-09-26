
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(pathname);
j=1;
Z = zeros(13,11);
for idx  = 3:length(directory)
    file = directory(idx).name;
    path = directory(idx).folder;
    A = readmatrix([path,'/',file]);
    Z(j,:) = A(end,:)';
    j = j+1;
end

for i=1:13
    dprime(i) = norminv(Z(i,3)/(Z(i,3)+Z(i,5))) - norminv(Z(i,6)/(Z(i,4)+Z(i,6)));
end

figure,plot(sort(dprime), 'Color',[0 0 0],'LineWidth',1.5);
xlabel('Sessions');
ylabel('Performance (d prime)');
set(gca,'TickDir','out','fontsize',14');

%%
nHits = 222; 
nCR = 96;
nMiss = 53;
nFA = 58;

dprime = norminv(nHits/(nHits+nMiss)) - norminv(nFA/(nFA+nCR));

