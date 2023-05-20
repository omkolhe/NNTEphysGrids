function [Encoder] = detectVelTrigStart(Encoder,velThreshold,tol,nSlope,nRejectWindow,LFPTimes)
 
disp('Thresholding for vel = '+ string(velThreshold) + ' with tolerance =' + string(tol));
Encoder.velTrig = find(Encoder.vel < velThreshold+tol & Encoder.vel >velThreshold-tol); % finding the points in time where velocity crosses the set threshold
padfordiff = nSlope; % number of points in front and back to calculate the instantaneous slope
for i=1:size(Encoder.velTrig,2) % Checking for slope at each points and rejecting negative slope points
    if (Encoder.velTrig(i)<padfordiff) || (Encoder.velTrig(i)>size(Encoder.vel,2)-padfordiff), Encoder.velTrig = NaN; end
    slope = mean(Encoder.vel(Encoder.velTrig(i):Encoder.velTrig(i)+padfordiff))-mean(Encoder.vel(Encoder.velTrig(i)-padfordiff:Encoder.velTrig(i)));
    if slope<0, Encoder.velTrig(i)=NaN; end
end
Encoder.velTrig = Encoder.velTrig(~isnan(Encoder.velTrig)); % rejecting all NaNs corresponding to negative slope
Encoder.goodInd = find(diff(Encoder.velTrig)>nRejectWindow); % rejecting points that are very close to each other (consequetive points within tolerance).
Encoder.velTrig = Encoder.velTrig(Encoder.goodInd);

Encoder.nTrig = numel(Encoder.velTrig);
for i=1:Encoder.nTrig
    Encoder.velTrig(2,i) = max(find(LFPTimes<=Encoder.time(Encoder.velTrig(1,i))));
end

disp('Number of velocity trigers found: '+ string(Encoder.nTrig))

end

