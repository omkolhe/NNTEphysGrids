hitLeverTraces = cell2mat(arrayfun(@(s) s.trace, IntanBehaviour.MIHitTrace,'UniformOutput',false));
FALeverTraces = cell2mat(arrayfun(@(s) s.trace, IntanBehaviour.MIFATrace,'UniformOutput',false));

figure('Name','Average Lever Traces for Cue Hits and FAs with allignment at MI');
hold on;
plot(IntanBehaviour.MIHitTrace(1).time,hitLeverTraces,'Color',[1 0 0 0.1],'LineWidth',1.5);
plot(IntanBehaviour.MIHitTrace(1).time,mean(hitLeverTraces,2),'Color',[1 0 0 1],'LineWidth',2);
yline(IntanBehaviour.MIcutoffHit,'--.b','MI Threshold','LabelHorizontalAlignment','left'); 

plot(IntanBehaviour.MIFATrace(1).time,FALeverTraces,'Color',[0 0 1 0.1],'LineWidth',1.5);
plot(IntanBehaviour.MIFATrace(1).time,mean(horzcat(IntanBehaviour.MIFATrace(1:end).trace),2),'Color',[0 0 1 1],'LineWidth',2);
yline(IntanBehaviour.MIcutoffFA,'--.b','MI Threshold','LabelHorizontalAlignment','left'); 
xline(0,'--r','MI','LabelVerticalAlignment','top');ylim([0 0.1]);
ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Hits and False Alarms');box off;

%% Truncating traces to just capture the rise

st = 1500-300;
sp = 1500+300;

hitLeverTraces = hitLeverTraces(st:sp,:);
FALeverTraces = FALeverTraces(st:sp,:);

figure('Name','Average Lever Traces for Cue Hits and FAs with allignment at MI');
hold on;
plot(st:1:sp,hitLeverTraces,'Color',[1 0 0 0.1],'LineWidth',1.5);
plot(st:1:sp,mean(hitLeverTraces,2),'Color',[1 0 0 1],'LineWidth',2);
yline(IntanBehaviour.MIcutoffHit,'--.b','MI Threshold','LabelHorizontalAlignment','left'); 

plot(st:1:sp,FALeverTraces,'Color',[0 0 1 0.1],'LineWidth',1.5);
plot(st:1:sp,mean(FALeverTraces,2),'Color',[0 0 1 1],'LineWidth',2);
yline(IntanBehaviour.MIcutoffFA,'--.b','MI Threshold','LabelHorizontalAlignment','left'); 
xline(0,'--r','MI','LabelVerticalAlignment','top');ylim([0 0.1]);
ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Hits and False Alarms');box off;


%% Performing PCA 

waveforms = [hitLeverTraces,FALeverTraces]';

figure();
plot(waveforms');
ylabel("Voltage (\mu V)")
xlabel("Time");

% covmatrix = corr(waveforms);
covmatrix = (waveforms'*waveforms);
covmatrix = covmatrix/size(waveforms,1);

figure();
imagesc(covmatrix);
colormap(jet);
colorbar;

[V,D] = eig(covmatrix);

q(:,1) = V(:,end);
q(:,2) = V(:,end-1);
q(:,3) = V(:,end-2);
q(:,4) = V(:,end-3);

figure();
plot(q);
ylabel("Voltage (\mu V)")
xlabel("Time");

projq = waveforms*q;

figure();
scatter3(projq(1:143,1),projq(1:143,2),projq(1:143,3),'b.','lineWidth',2);
hold on;
scatter3(projq(144:end,1),projq(144:end,2),projq(144:end,3),'r.','lineWidth',2);
ylabel("bk")
xlabel("ak");

thresA = 0;

for i=1:size(spikes,2)
    if projq(i,1) <= thresA %&& projq(i,2)< 0
        projq(i,3) = 1;
    elseif projq(i,1) > thresA %&& projq(i,2)< 75
        projq(i,3) = 2;
    else 
        projq(i,3) = 3;
    end 
end

figure(),hold on;
for i=1:size(spikes,2)
    if projq(i,3) == 1
        plot(projq(i,1),projq(i,2),'b.','lineWidth',2);
    end
    if projq(i,3) == 2
        plot(projq(i,1),projq(i,2),'g.','lineWidth',2);
    end
    if projq(i,3) == 3
        plot(projq(i,1),projq(i,2),'r.','lineWidth',2);
    end
end


figure();
subplot(3,1,1)
for i=1:size(spikes,2)
    if projq(i,3) == 1
        plot(waveforms(i,:),'Color',[0 0 1 0.1],'lineWidth',1);
        hold on;
    end  
end 
subplot(3,1,2)
for i=1:size(spikes,2)
    if projq(i,3) == 2
        plot(waveforms(i,:),'Color',[0 1 0 0.1],'lineWidth',1);
        hold on;
    end  
end
subplot(3,1,3)
for i=1:size(spikes,2)
    if projq(i,3) == 3
        plot(waveforms(i,:),'Color',[1 0 0 0.1],'lineWidth',1);
        hold on;
    end  
end
ylabel("Voltage (\mu V)");
xlabel("Time");

ak = projq;
ak(:,2) = [];
bk = projq;
bk(:,1) = [];

a2 = mean(ak(ak(:,2)==2)); 
a1 = mean(ak(ak(:,2)==1)); 

b2 = mean(bk(bk(:,2)==2)); 
b1 = mean(bk(bk(:,2)==1)); 

spike2 = a2*q(:,1) + b2*q(:,2);
spike1 = a1*q(:,1) + b1*q(:,2);

figure();
plot(spike1);
hold on;
 plot(spike1,'r','lineWidth',2);

 figure();
plot(spike2);
hold on;
 plot(spike2,'r','lineWidth',2);


