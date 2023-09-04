function [PSDCh, f] = getAvgPSD(behaviourTrace,parameters)

for trialno = 1:size(behaviourTrace,2)
    xf1 = behaviourTrace(trialno).rawLFP;
    for i=1:parameters.rows
        for j=1:parameters.cols
            [PSDCh(trialno,(i-1)*parameters.cols + j,:) ,f] = pwelch(squeeze(xf1(parameters.rows,parameters.cols,:)),600,0,1000,parameters.Fs);
        end
    end
end

