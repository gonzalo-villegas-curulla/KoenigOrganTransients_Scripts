clc, clear;
files = dir('A*.mat');

for IDX = 1 : length(files)
load(files(IDX).name);
x = PR_data.pressureData;
for idx = 1:min(size(x))
    x(idx,:) = x(idx,:) - mean(x(idx,1:1e3),'omitnan');
end


figure(1); clf;
plot(x(1,:));
hold on;
plot(x(2,:));
plot(x(3,:));
plot(x(4,:));
plot(x(5,:));
yyaxis right;
plot(x(6,:));
yyaxis left;
grid on;
title(sprintf(files(IDX).name));

pause();

end

%%
dt = 1/51.2e3;
len = length(x);
figure(10); clf; hold on;
for idx = 1 : 5
    plot( [0:len-1]*dt, x(idx,:) );
end
xlabel('Time [s]');
ylabel('Pressure [Pa]');
grid on; box on;

yyaxis right;
plot([0:len-1]*dt, x(6,:));
ylabel('Key velocity [m/s]');

%%
