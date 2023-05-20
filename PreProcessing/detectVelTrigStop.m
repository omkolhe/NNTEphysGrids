function [Encoder] = detectVelTrigStop(Encoder,velThreshold,tol,nSlope,nRejectWindow,LFPTimes)
 
disp('Thresholding for vel = '+ string(velThreshold) + ' with tolerance =' + string(tol));
Encoder.velTrigStop = find(Encoder.vel < velThreshold+tol & Encoder.vel >velThreshold-tol); % finding the points in time where velocity crosses the set threshold
padfordiff = nSlope; % number of points in front and back to calculate the instantaneous slope
for i=1:size(Encoder.velTrigStop,2) % Checking for slope at each points and rejecting negative slope points
    if (Encoder.velTrigStop(i)<padfordiff) || (Encoder.velTrigStop(i)>size(Encoder.vel,2)-padfordiff), Encoder.velTrigStop = NaN; end
    slope = mean(Encoder.vel(Encoder.velTrigStop(i):Encoder.velTrigStop(i)+padfordiff))-mean(Encoder.vel(Encoder.velTrigStop(i)-padfordiff:Encoder.velTrigStop(i)));
    if slope>0, Encoder.velTrigStop(i)=NaN; end
end
Encoder.velTrigStop = Encoder.velTrigStop(~isnan(Encoder.velTrigStop)); % rejecting all NaNs corresponding to negative slope
Encoder.goodIndStop = find(diff(Encoder.velTrigStop)>nRejectWindow); % rejecting points that are very close to each other (consequetive points within tolerance).
Encoder.velTrigStop = Encoder.velTrigStop(Encoder.goodIndStop);

Encoder.nTrigStop = numel(Encoder.velTrigStop);
for i=1:Encoder.nTrigStop
    Encoder.velTrigStop(2,i) = max(find(LFPTimes<=Encoder.time(Encoder.velTrigStop(1,i))));
end

disp('Number of velocity trigers found: '+ string(Encoder.nTrigStop))

end

