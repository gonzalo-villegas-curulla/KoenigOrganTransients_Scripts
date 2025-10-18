load('data_proc.mat');

pipelist = [1:22];

rho = 1.2; co = 340; co2 = co*co; 

Ta_geo = []; Tb_geo = []; Tc_geo = []; Td_geo = [];
for idx = 1 :length(pipelist)

    sample_select = pipelist(idx);

    P0 = data_proc.Ppall(sample_select);
    
    Vg = data_proc.Vgrv(sample_select);
    Vf = data_proc.Vf(sample_select);

    Spall_geo = data_proc.Spall_geom(sample_select);
    Sin_geo   = data_proc.Sin_geom(sample_select);
    Sj_geo    = data_proc.Sjet_geom(sample_select);


    Ta_geo = [Ta_geo, Vg/(Spall_geo*co2)*sqrt(P0/rho)];
    Tb_geo = [Tb_geo, Vg/(Sin_geo*co2)*sqrt(P0/rho)];
    Tc_geo = [Tc_geo, Vf/(Sin_geo*co2)*sqrt(P0/rho)];
    Td_geo = [Td_geo, Vf/(Sj_geo*co2)*sqrt(P0/rho)];   

end

fax = 12*log2(data_proc.F1/440);
LW = 1.5;

figure(1); clf; hold on; grid on; box on;
%
plot(fax, 1e3*Ta_geo, 'r-s', 'DisplayName','T_a^{geo}','linewidth',LW);
plot(fax, 1e3*data_proc.Amax, 'r-d', 'DisplayName','A_{eff}','linewidth',LW,...
    'markerfacecolor','r');
%
plot(fax, 1e3*Tb_geo, 'b-s', 'DisplayName','T_b^{geo}','linewidth',LW);
plot(fax, 1e3*data_proc.B, 'b-d', 'DisplayName','B_{eff}','linewidth',LW,...
    'markerfacecolor','b');
%
plot(fax, 1e3*Tc_geo, 'g-s', 'DisplayName','T_c^{geo}','linewidth',LW);
plot(fax, 1e3*data_proc.C, 'g-d', 'DisplayName','C_{eff}','linewidth',LW,...
    'markerfacecolor','g');
%
plot(fax, 1e3*Td_geo, 'm-s', 'DisplayName','T_d^{geo}','linewidth',LW);
plot(fax, 1e3*data_proc.D, 'm-d', 'DisplayName','D_{eff}','linewidth',LW,...
    'markerfacecolor','m');
%
xlabel('12log_2(F1/440)');
ylabel('ms');
legend('show', 'location','northwest');