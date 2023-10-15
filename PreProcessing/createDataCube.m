function [LFP] = createDataCube(LFP,rows,cols,chMap)
    
LFP.LFPdatacube = NaN(rows,cols,size(LFP.LFP,2));

for i=1:rows
    for j=1:cols
        if(~ismember((cols*(i-1)+j),chMap))
            LFP.LFPdatacube(i,j,:) = LFP.LFP(cols*(i-1)+j,:);
        end
    end
end

end

