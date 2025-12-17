clc, clear;
load data_proc;

P0 = data_proc.Ppall_mean(:);
rho = 1.2;
co = 340; co2 = co*co;
% ============================================
Vg = data_proc.Vgrv(:); % Measured geometry
Vf = data_proc.Vf(:);   % Measured geometry
Sp_geo  = data_proc.Spall_geom(:); % Measured
Sin_geo = data_proc.Sin_geom(:); % Measured
Sj_geo  = data_proc.Sjet_geom(:); % Neasured
% ============================================
% (1)
SIG = 0;
% (2)
PRTg = data_proc.PRTgrv_mean(:); % Measured/deduced time
Ta = PRTg./1.45; % Model

% (3) 
Pg_targ = data_proc.Pgrv_mean(:); % Measured pressure
Pf_targ = data_proc.Pf_mean(:);   % Measured pressure

    ABchunk = (1-Pg_targ(:)./P0(:))./(Pg_targ(:)./P0(:).*(1-SIG^2) - Pf_targ./P0(:) + SIG^2);
    CDchunk = sqrt(Pf_targ./P0(:) ./ (Pg_targ./P0(:)*(1-SIG^2)-Pf_targ./P0(:) + SIG^2));

Tb = Ta(:)./sqrt(ABchunk);
Tc = Tb.*Vf./Vg; 
Td = Tc.*CDchunk;
% ============================================
Speff    = Vg(:)./Ta./co2.*sqrt(P0/rho);
Sin_eff  = Vg(:)./Tb./co2.*sqrt(P0/rho);
Sj_eff   = Vf(:)./Td./co2.*sqrt(P0/rho);

Sp_ratio  = Speff(:)./Sp_geo(:);
Sin_ratio = Sin_eff(:)./Sin_geo(:);
Sj_ratio  = Sj_eff(:)./Sj_geo(:);
% ============================================
figure(1); clf;
fax = 12*log2(data_proc.F1(:)/440);
plot(fax, (Sp_ratio), '-sk', 'markeredgecolor','k');
hold on;
plot(fax, (Sin_ratio), '-ok', 'markeredgecolor','k','markerfacecolor','k');
plot(fax, (Sj_ratio), '-^k','markeredgecolor','k','markerfacecolor','k');
grid on;

xlabel('$12 \times log_2(f_1/440 \ Hz)$', 'interpreter','latex');
ylabel('$S/S^{geo}$ [n.u.]','interpreter','latex');
ax=gca;
ax.YScale = 'log';
legend('$S_p^{ratio}, \ \Sigma = 0$',...
       '$S_{in}^{ratio}, \ \Sigma = 0$',...
       '$S_j^{ratio}, \ \Sigma = 0$', 'location','best','interpreter','latex');
