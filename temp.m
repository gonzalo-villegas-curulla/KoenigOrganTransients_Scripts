
load('A07_pitchFsharp_20240327_181408.mat');
x = PR_data.pressureData;
fs=51.2e3; dt=1/fs;
for idx = 1:6
    x(idx,:) = x(idx,:) - mean(x(idx, 1:1e3));
end


ll = 100e3 : 3e6;
frag = x(1, ll);
flims = [3, 100];

figure(1); clf;
plot([0:length(frag)-1]*dt,frag );
ax(1) = gca;

figure(2); clf;
cwt(frag, fs, 'FrequencyLimits', flims,...
    'VoicesPerOctave', 30);
caxis([1 6]);
ax(2) = gca;

linkaxes(ax, 'x');
