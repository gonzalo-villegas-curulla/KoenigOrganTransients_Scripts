files = dir('A*.mat');
load(files(2).name);
x = PR_data.pressureData;
for idx = 1 : 6
    x(idx,:) = x(idx,: ) - mean(x(idx,1:1e3));
end


%%
pres = x(2,:);

fs = 51.2e3; dt = 1/fs;
t = [0:length(pres)-1]*dt;

figure; plot(t,pres);

fac = 10;
fs = fs/fac;

pres = resample(pres, fs, fac*fs);
t=resample(t, fs, 10*fs);
hold on; plot(t,pres);


%%
seg = pres(389151:509137);
seg = pres(389151:447228);

fois = 1e-2:0.1:80;
sord = [1,30];

Cdata = aslt(seg, fs, fois, 19, srord, 0);
% Cdata = aslt(xSignal, fs, fois, 5, srord, 0);


figure(10); clf; 
imagesc(t, fois, log(Cdata));
set(gca, 'ydir','normal');
colormap jet; colorbar;title('power (log)');
xlabel('Time [frames]');ylabel('Freq [Hz]');
caxis([-10 10])


figure(11); clf; 
imagesc(t, fois, Cdata);
set(gca, 'ydir','normal');
xlabel('Time [frames]');ylabel('Freq [Hz]');
title('power (lin)');
colormap jet; colorbar;
% caxis([1e-2 10]);
caxis([1e-2 20])

%%


%%