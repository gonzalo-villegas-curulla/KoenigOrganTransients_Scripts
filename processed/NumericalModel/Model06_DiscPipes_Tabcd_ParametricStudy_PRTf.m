clc, clear;
load data_proc.mat

% %%%%%%%%%%%%%%%%   CUSTOM USER PARAMETERS  %%%%%%%%%%%%%%%%%%%

pipelist = [8, 18]; % Sampled pipes [1:22]= ["A01",..., "A56"]inventory

% Pallet valve openig time 

ValveRampInit = 0.100;  % [s] T-Start opening-ramp pallet valve
ValveRampEnd  = 0.101;  % [s] T-Finish opening-ramp (DROPIC robot time)

NUMVALS = 10;

% %%%%%%%%%%%%%%   END OF CUSTOM USER PARAMETERS  %%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------------------
%  1   2   3   4   5   6   7   8   9   10   11   12   13   14   15   16     <=== PipeNum
%          1   2   3   4   5       6    7    8         9        10          <=== Sample Num
%          X   X   X   X   X       X    X    X         X         X          <=== Trusted
% ------------------------------------------------------------------------
%  17   18   19   20   21   22   23   24   25   26   27   28   29   30
%  11        12                       13   14        15        16
%   X         X                        X    X
% ------------------------------------------------------------------------
%  31   32   33   34   35   36   37   38   39   40   41   42   43   44 ...
%       17        18             19        20        21             22
%                  X
% ------------------------------------------------------------------------



% ========= Physical constants =======
rho = 1.2;
co  = 340;          co2 = co^2;
P0  = 820;

% ======= Simulation parameters ==============
fs   = 51.2e3; 
dt   = 1/fs;
Tend = 0.300;
tvec = [0:dt:Tend]';


sample_select = pipelist(:);

% PIPE 1 refs

Ta_ref1 = data_proc.Amax(sample_select(1));
Tb_ref1 = data_proc.B(sample_select(1));
Tc_ref1 = data_proc.C(sample_select(1));
Td_ref1 = data_proc.D(sample_select(1));
Tomega_ref = 1e-3;

% Pipe 2 refs
Ta_ref2 = data_proc.Amax(sample_select(2));
Tb_ref2 = data_proc.B(sample_select(2));
Tc_ref2 = data_proc.C(sample_select(2));
Td_ref2 = data_proc.D(sample_select(2));

% Parametric slow variation range
Ta_log     = 1e-3*logspace(-0.15,1.31,NUMVALS);
Tb_log     = 1e-3*logspace(0.31, 1.7, NUMVALS);
Tc_lin     = 1e-3*linspace(1.0,2.0, NUMVALS);
Td_lin     = 1e-3*linspace(1.0, 2.0, NUMVALS);
Tomega_log = 1e-3*logspace(-1.3, 1.31, 10);

   


    % ===================================
    %           Plot results
    % ===================================

% FUNCTION: [PRTg] = parametrify(Ta, Tb, Tc, Td, Tomega, SIGMA, ValveRampInit, Tend, tvec)    

% ============= Ta ===============

            SIGl = 0.0;
            SIGh = 0.5;

fig1 = figure(1); clf; ax1 = axes(fig1); hold on; box on; grid on;
% PIPE 1
resultL = parametrifyMULT_func(Ta_log, Tb_ref1, Tc_ref1,Td_ref1, Tomega_ref, SIGl, ValveRampInit, Tend, tvec);
resultH = parametrifyMULT_func(Ta_log, Tb_ref1, Tc_ref1,Td_ref1, Tomega_ref, SIGh, ValveRampInit, Tend, tvec);
plot(ax1, 1e3*Ta_log(find(resultL)), 1e3*resultL, '-o');
plot(ax1, 1e3*Ta_log(find(resultH)), 1e3*resultH, '-o');
%
% PIPE 2
resultL = parametrifyMULT_func(Ta_log, Tb_ref2, Tc_ref2,Td_ref2, Tomega_ref, SIGl, ValveRampInit, Tend, tvec);
resultH = parametrifyMULT_func(Ta_log, Tb_ref2, Tc_ref2,Td_ref2, Tomega_ref, SIGh, ValveRampInit, Tend, tvec);
plot(ax1, 1e3*Ta_log(find(resultL)), 1e3*resultL, '-o');
plot(ax1, 1e3*Ta_log(find(resultH)), 1e3*resultH, '-o');
%
xlabel('$\mathcal{T}_a$ [ms]','interpreter','latex');
xlim([0.5 21]);
ylabel('PRT$_f$ [ms]', 'interpreter','latex');
legend({'$A_3^{\#}, \ \Sigma_L$','$A_3^{\#}, \ \Sigma_H$','$A_5^{\#}, \ \Sigma_L $','$A_5^{\#}, \ \Sigma_H$'},...
    'interpreter','latex', 'location','northwest');
title(sprintf('SIGMA_{L,H}= %1.1f, %1.1f', SIGl, SIGh));
ax=gca; ax.XLim(1) = 0;
ax.YLim(1) = 0;
% ============= Tb ===============

            SIGl = 0.0;
            SIGm = 0.5;
            SIGh = 0.9;

fig2 = figure(2); clf; ax2 = axes(fig2); hold on; box on; grid on;
% PIPE 1
resultL = parametrifyMULT_func(Ta_ref1, Tb_log, Tc_ref1,Td_ref1, Tomega_ref, SIGl, ValveRampInit, Tend, tvec);
resultH = parametrifyMULT_func(Ta_ref1, Tb_log, Tc_ref1,Td_ref1, Tomega_ref, SIGh, ValveRampInit, Tend, tvec);
plot(ax2, 1e3*Tb_log(find(resultL)), 1e3*resultL, '-o');
plot(ax2, 1e3*Tb_log(find(resultH)), 1e3*resultH, '-o');
%
% PIPE 2
resultL = parametrifyMULT_func(Ta_ref2, Tb_log, Tc_ref2,Td_ref2, Tomega_ref, SIGl, ValveRampInit, Tend, tvec);
resultM = parametrifyMULT_func(Ta_ref2, Tb_log, Tc_ref2,Td_ref2, Tomega_ref, SIGm, ValveRampInit, Tend, tvec);
resultH = parametrifyMULT_func(Ta_ref2, Tb_log, Tc_ref2,Td_ref2, Tomega_ref, SIGh, ValveRampInit, Tend, tvec);
plot(ax2, 1e3*Tb_log(find(resultL)), 1e3*resultL, '-o');
plot(ax2, 1e3*Tb_log(find(resultM)), 1e3*resultM, '-o');
plot(ax2, 1e3*Tb_log(find(resultH)), 1e3*resultH, '-o');
%
xlabel('$\mathcal{T}_b$ [ms]','interpreter','latex');
xlim([1 51]);
ylabel('PRT$_f$ [ms]', 'interpreter','latex');
legend({'$A_3^{\#}, \ \Sigma_L$','$A_3^{\#}, \ \Sigma_H$','$A_5^{\#}, \ \Sigma_L $','$A_5^{\#}, \ \Sigma_M $','$A_5^{\#}, \ \Sigma_H$'},...
    'interpreter','latex', 'location','best');
title(sprintf('SIGMA_{L,M,H}= %1.1f, %1.1f ,%1.1f', SIGl, SIGm, SIGh));
ax=gca; ax.XLim(1) = 0;

% ============= Tc ===============

            SIGl = 0.0;
            SIGm = 0.5;
            SIGh = 0.9;

fig3 = figure(3); clf; ax3 = axes(fig3); hold on; box on; grid on;
% PIPE 1
resultL = parametrifyMULT_func(Ta_ref1, Tb_ref1, Tc_lin,Td_ref1, Tomega_ref, SIGl, ValveRampInit, Tend, tvec);
resultH = parametrifyMULT_func(Ta_ref1, Tb_ref1, Tc_lin,Td_ref1, Tomega_ref, SIGh, ValveRampInit, Tend, tvec);
plot(ax3, 1e3*Tc_lin(find(resultL)), 1e3*resultL, '-o');
plot(ax3, 1e3*Tc_lin(find(resultH)), 1e3*resultH, '-o');
%
% PIPE 2
resultL = parametrifyMULT_func(Ta_ref2, Tb_ref2, Tc_lin,Td_ref2, Tomega_ref, SIGl, ValveRampInit, Tend, tvec);
resultM = parametrifyMULT_func(Ta_ref2, Tb_ref2, Tc_lin,Td_ref2, Tomega_ref, SIGm, ValveRampInit, Tend, tvec);
resultH = parametrifyMULT_func(Ta_ref2, Tb_ref2, Tc_lin,Td_ref2, Tomega_ref, SIGh, ValveRampInit, Tend, tvec);
plot(ax3, 1e3*Tc_lin(find(resultL)), 1e3*resultL, '-o');
plot(ax3, 1e3*Tc_lin(find(resultM)), 1e3*resultM, '-o');
plot(ax3, 1e3*Tc_lin(find(resultH)), 1e3*resultH, '-o');
%
xlabel('$\mathcal{T}_c$ [ms]','interpreter','latex');
xlim([0.9 2.1]);
ylabel('PRT$_f$ [ms]', 'interpreter','latex');
legend({'$A_3^{\#}, \ \Sigma_L$','$A_3^{\#}, \ \Sigma_H$','$A_5^{\#}, \ \Sigma_L $','$A_5^{\#}, \ \Sigma_M $','$A_5^{\#}, \ \Sigma_H$'},...
    'interpreter','latex', 'location','best');
title(sprintf('SIGMA_{L,M,H}= %1.1f, %1.1f ,%1.1f', SIGl, SIGm, SIGh));
ax=gca; ax.YLim(1) = 0;

% ============= Td ===============

            SIGl = 0.0;
            SIGm = 0.5;
            SIGh = 0.9;

fig4 = figure(4); clf; ax4 = axes(fig4); hold on; box on; grid on;
% PIPE 1
resultL = parametrifyMULT_func(Ta_ref1, Tb_ref1, Tc_ref1,Td_lin, Tomega_ref, SIGl, ValveRampInit, Tend, tvec);
resultH = parametrifyMULT_func(Ta_ref1, Tb_ref1, Tc_ref1,Td_lin, Tomega_ref, SIGh, ValveRampInit, Tend, tvec);
plot(ax4, 1e3*Td_lin(find(resultL)), 1e3*resultL, '-o');
plot(ax4, 1e3*Td_lin(find(resultH)), 1e3*resultH, '-o');
%
% PIPE 2
resultL = parametrifyMULT_func(Ta_ref2, Tb_ref2, Tc_ref2,Td_lin, Tomega_ref, SIGl, ValveRampInit, Tend, tvec);
resultM = parametrifyMULT_func(Ta_ref2, Tb_ref2, Tc_ref2,Td_lin, Tomega_ref, SIGm, ValveRampInit, Tend, tvec);
resultH = parametrifyMULT_func(Ta_ref2, Tb_ref2, Tc_ref2,Td_lin, Tomega_ref, SIGh, ValveRampInit, Tend, tvec);
plot(ax4, 1e3*Td_lin(find(resultL)), 1e3*resultL, '-o');
plot(ax4, 1e3*Td_lin(find(resultM)), 1e3*resultM, '-o');
plot(ax4, 1e3*Td_lin(find(resultH)), 1e3*resultH, '-o');
%
xlabel('$\mathcal{T}_d$ [ms]','interpreter','latex');
xlim([0.9 2.1]);
ylabel('PRT$_f$ [ms]', 'interpreter','latex');
legend({'$A_3^{\#}, \ \Sigma_L$','$A_3^{\#}, \ \Sigma_H$','$A_5^{\#}, \ \Sigma_L $','$A_5^{\#}, \ \Sigma_M $','$A_5^{\#}, \ \Sigma_H$'},...
    'interpreter','latex', 'location','best');
title(sprintf('SIGMA_{L,M,H}= %1.1f, %1.1f ,%1.1f', SIGl, SIGm, SIGh));
ax=gca; ax.YLim(1) = 0;


% ============= Tau Omega ===============

            SIGl = 0.0;
            SIGm = 0.5;
            SIGh = 0.9;

fig5 = figure(5); clf; ax5 = axes(fig5); hold on; box on; grid on;
% PIPE 1
resultL = parametrifyMULT_func(Ta_ref1, Tb_ref1, Tc_ref1,Td_ref1, Tomega_log, SIGl, ValveRampInit, Tend, tvec);
resultH = parametrifyMULT_func(Ta_ref1, Tb_ref1, Tc_ref1,Td_ref1, Tomega_log, SIGh, ValveRampInit, Tend, tvec);
plot(ax5, 1e3*Tomega_log(find(resultL)), 1e3*resultL, '-o');
plot(ax5, 1e3*Tomega_log(find(resultH)), 1e3*resultH, '-o');
%
% PIPE 2
resultL = parametrifyMULT_func(Ta_ref2, Tb_ref2, Tc_ref2,Td_ref2, Tomega_log, SIGl, ValveRampInit, Tend, tvec);
resultM = parametrifyMULT_func(Ta_ref2, Tb_ref2, Tc_ref2,Td_ref2, Tomega_log, SIGm, ValveRampInit, Tend, tvec);
resultH = parametrifyMULT_func(Ta_ref2, Tb_ref2, Tc_ref2,Td_ref2, Tomega_log, SIGh, ValveRampInit, Tend, tvec);
plot(ax5, 1e3*Tomega_log(find(resultL)), 1e3*resultL, '-o');
plot(ax5, 1e3*Tomega_log(find(resultM)), 1e3*resultM, '-o');
plot(ax5, 1e3*Tomega_log(find(resultH)), 1e3*resultH, '-o');
%
xlabel('$\tau_{\Omega}$ [ms]','interpreter','latex');
xlim([0.01 21]);
ylabel('PRT$_f$ [ms]', 'interpreter','latex');
legend({'$A_3^{\#}, \ \Sigma_L$','$A_3^{\#}, \ \Sigma_H$','$A_5^{\#}, \ \Sigma_L $','$A_5^{\#}, \ \Sigma_M $','$A_5^{\#}, \ \Sigma_H$'},...
    'interpreter','latex', 'location','best');
title(sprintf('SIGMA_{L,M,H}= %1.1f, %1.1f ,%1.1f', SIGl, SIGm, SIGh));


% obj=findobj('type','line'); for idx=1:length(obj) obj(idx).LineWidth=2; end

