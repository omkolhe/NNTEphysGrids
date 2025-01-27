function [Z,goodChMap,badCh] = readImp(imppath,impfile,electrode_map,cutoff)


A = readmatrix([imppath,'/',impfile]);
Z = A(:,5);
Z = Z(electrode_map,1);

badCh = find(Z>cutoff);

chMap = linspace(1,32,32);
[X,Y] = ismember(badCh,chMap);
chMap(Y(X)) = [];
goodChMap = chMap';

disp(['Number of bad channels based of impedance : ', num2str(numel(badCh))]);

end

