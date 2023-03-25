function [Intan] = removeBadCh(Intan,badCh)
%REMOVEBADCH Summary of this function goes here
%   Detailed explanation goes here
    chMap = linspace(1,32,32);
    [X,Y] = ismember(badCh,chMap);
    chMap(Y(X)) = [];
    Intan.goodChMap = chMap';
    Intan.badChMap = badCh;
    Intan.allIntanBad = Intan.allIntan(badCh,:);
    Intan.allIntan(badCh,:) = [];
end

