clc, clear;

load data_proc.mat

pipelist = [1:22]; % \in[1,22] (DO NOT INCLUDE MORE THAN 3 PIPES at a time)

% Pallet valve openig time 
ValveRampInit = 0.100; % [s] T-Start opening-ramp pallet valve
ValveRampEnd  = 0.10001;  % [s] T-Finish opening-ramp (DROPIC robot time)


% ========= Physical constants =======
rho = 1.2;
co  = 340;          co2 = co^2;

% ======= Simulation parameters ==============
fs   = 2*51.2e3; 
dt   = 1/fs;
Tend = 0.300;
tvec = [0:dt:Tend]';


% Parameter value allocation and simulation run:
MX_results = zeros(length(pipelist), 5); % Ppall, Pgrv, Pf, PRTgrv, PRTf

res_Pgrv_trg = zeros(length(pipelist), 1);
res_Pf_trg   = zeros(length(pipelist), 1);
res_PRTgrv   = zeros(length(pipelist), 1);
res_PRTf     = zeros(length(pipelist), 1);
Ares = [];
Bres = [];

% Perform simulation ==================================
for PLidx = 1:length(pipelist) 

    sample_select = pipelist(PLidx);

    Xi     = data_proc.Amax(sample_select)/data_proc.B(sample_select);
    Zeta   = data_proc.C(sample_select)/data_proc.D(sample_select);
    Volrat = data_proc.B(sample_select)/data_proc.C(sample_select);

    Atmp = 0.8*data_proc.PRTgrv_mean(sample_select); % From PRT   
    Btmp = Atmp/Xi;    % => data_proc.B(sample_select)./data_proc.Amax(sample_select)*Atmp;
    Ctmp = Btmp/Volrat;% => data_proc.C(sample_select)./data_proc.B(sample_select)*Btmp;
    Dtmp = Ctmp/Zeta;  % => data_proc.D(sample_select)./data_proc.C(sample_select)*Ctmp;

    
    Spall_eff(PLidx)  = data_proc.Vgrv(PLidx)/(Atmp*co2)*sqrt(data_proc.Ppall(PLidx)/rho);
    Spall_geom(PLidx) = data_proc.Spall_geom(PLidx);

    Sin_eff(PLidx)  = data_proc.Vgrv(PLidx)/(Btmp*co2)*sqrt(data_proc.Ppall(PLidx)/rho);
    Sin_eff2(PLidx) = data_proc.Vf(PLidx)/(Ctmp*co2)*sqrt(data_proc.Ppall(PLidx)/rho); 
    Sin_geom(PLidx) = data_proc.Sin_geom(PLidx);
    
    Sj_geom(PLidx) = data_proc.Sjet_geom(PLidx);
    Sj_eff(PLidx)  = data_proc.Vf(PLidx)/(Dtmp*co2)*sqrt(data_proc.Ppall(PLidx)/rho);
    % % % Sj_eff(PLidx)  = Sj_geom(PLidx);
    
end


    % ===================================
    %           Plot results
    % ===================================

figure(12); clf;

plot(1e6*Spall_geom, 1e6*Spall_eff, '-o');
xlabel('Spall geom [mm^2]');ylabel('Spall eff [mm^2]'); grid on;
xlim([0 1200]);
hold on; plot([0 600],[0 600],'-k'); 
axis equal; title('Spall');
ax=gca;ax.XLim(1) = 0; ax.YLim(1) = 0;
%
figure(13); clf;
plot(1e6*Sin_geom, 1e6*Sin_eff, '-o');
hold on;grid on; 
plot(1e6*Sin_geom, 1e6*Sin_eff2, '-o');
xlabel('Sin geom [mm^2]'); ylabel('Sin eff [mm^2]');
plot([0 50],[0 50],'-k');
axis equal; ylim([0 50]);
title('Sin')
legend('Vgrv and B','Vf and C', 'location', 'best');
ax=gca; ax.XLim(1) = 0; ax.YLim(1) = 0;
%
figure(14); clf;
plot(1e6*Sj_geom, 1e6*Sj_eff, '-o');
hold on;
plot([0 50],[0 50],'-k');
grid on; title('Sjet');
axis equal;
ax=gca;
ax.YLim(1) = 0;
ax.XLim(1) = 0;
xlabel('Sjet geom [mm^2]');ylabel('Sjet eff [mm^2]');

