clc; clear;
load data_proc.mat;
lp  = 5e-3;
lin = 1e-3;
lj  = 1e-3;

rho = 1.2;
co = 340; co2 = co*co;

P0 = data_proc.Ppall';

Vg = data_proc.Vgrv';
Vf = data_proc.Vf';

Spinf = data_proc.Spall_geom';
Sin   = data_proc.Sin_geom';
Sj    = data_proc.Sjet_geom';



tp = lp.*sqrt(rho./P0);
tin = lin.*sqrt(rho./P0);
tj = lj.*sqrt(rho./P0);

Ta = Vg./(Spinf*co2).*sqrt(P0/rho);
Tb = Vg./(Sin*co2).*sqrt(P0/rho);
Tc = Vf./(Sin*co2).*sqrt(P0/rho);
Td = Vf./(Sj*co2).*sqrt(P0/rho);


fax = 12*log2(data_proc.F1'/440);
figure(1); clf;
hold on; grid on; box on;
xlabel('12log_2(F_1/440)');
ylabel('ch.t. [ms]');
plot(fax, 1e3*tp, '-o');
plot(fax, 1e3*tin, '-o');
plot(fax, 1e3*tj, '-o');
plot(fax, 1e3*Ta, '-o');
plot(fax, 1e3*Tb, '-o');
plot(fax, 1e3*Tc, '-o');
plot(fax, 1e3*Td, '-o');
%
title('geometrical');
str = {'$\tau_p$', '$\tau_{in}$', '$\tau_j$', '$\mathcal{T}_a$',...
    '$\mathcal{T}_b$', '$\mathcal{T}_c$', '$\mathcal{T}_d$'};
legend(str, 'location','best', 'interpreter','latex', 'fontsize', 14);
ax=gca; ax.YScale = 'log';
lax = findobj('type','line');
for idx = 1:length(lax)
    lax(idx).LineWidth = 2;
end


%

figure(2); clf; hold on; box on; grid on;
plot(fax, 1e3*tp, '-o');
plot(fax, 1e3*tin, '-o');
plot(fax, 1e3*tj, '-o');
plot(fax, 1e3*data_proc.Amax', '-o');
plot(fax, 1e3*data_proc.B', '-o');
plot(fax, 1e3*data_proc.C', '-o');
plot(fax, 1e3*data_proc.D', '-o');
title('effective values');
legend(str, 'location','best', 'interpreter','latex', 'fontsize', 14);
ax=gca; ax.YScale = 'log';
lax = findobj('type','line');
for idx = 1:length(lax)
    lax(idx).LineWidth = 2;
end
ylabel('ms');
xlabel('12log_2(F_1/440)');