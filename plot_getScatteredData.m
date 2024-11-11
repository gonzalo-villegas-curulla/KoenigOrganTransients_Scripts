ax = gca;
N = length(ax.Children);

X = [1,1];
Y = [1,1];
for idx = 1 :N
    X = [X, ax.Children(idx).XData];
    Y = [Y, ax.Children(idx).YData];
end

X = X(3:end);
Y = Y(3:end);

figure();
scatter(log(X),log(Y));


ax=gca;ax.XScale='log';
ax.YScale='log';
% xlim([8e-1,20]);
% ylim(1e3*[0.5063,5.4008]);
grid on;