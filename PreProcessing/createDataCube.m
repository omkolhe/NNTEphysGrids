function [LFP] = createDataCube(LFP,rows,cols,chMap)
    
LFP.LFPdatacube = NaN(rows,cols,size(LFP.LFP,2));

for i=1:rows
    for j=1:cols
        if(ismember((6*(i-1)+j),ChMap))
            LFP.LFPdatacube(i,j,:) = LFP.LFP(6*(i-1)+j,:);
        end
    end
end

end

