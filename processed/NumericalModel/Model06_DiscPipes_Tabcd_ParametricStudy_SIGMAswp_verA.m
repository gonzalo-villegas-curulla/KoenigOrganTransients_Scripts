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


% ============= Tb ===============


SIGvec = [0.0 : 0.05 : 1.0]';

fig2 = figure(2); clf; ax2 = axes(fig2); hold on; box on; grid on;
Tb_min     = min(data_proc.B);
Tb_sample8 = data_proc.B(8);

for SIG = [0.0 : 0.05 : 1.0]

    % PIPE 1 (def. 8)
    [PRTg1,~] = parametrifyMULT_func(Ta_ref1, Tb_log, Tc_ref1,Td_ref1, Tomega_ref, SIG, ValveRampInit, Tend, tvec);
    
    try
    plot(ax2, 1e3*Tb_log, 1e3*PRTg1, 'ok');
    end    

    
    % PIPE 2 (def. 18)
    
    [PRTg2,~] = parametrifyMULT_func(Ta_ref2, Tb_log, Tc_ref2,Td_ref2, Tomega_ref, SIG, ValveRampInit, Tend, tvec);

    try
    plot(ax2, 1e3*Tb_log, 1e3*PRTg2, 'dr');
    end

    drawnow();

end
%



xlabel('$\mathcal{T}_b$ [ms]','interpreter','latex');
xlim([1 51]);
ylabel('PRT$_g$ [ms]', 'interpreter','latex');
legend({'$A_3^{\#}, \ \Sigma_L$','$A_3^{\#}, \ \Sigma_H$','$A_5^{\#}, \ \Sigma_L $','$A_5^{\#}, \ \Sigma_M $','$A_5^{\#}, \ \Sigma_H$'},...
    'interpreter','latex', 'location','best');
% title(sprintf('SIGMA_{L,M,H}= %1.1f, %1.1f ,%1.1f', SIGl, SIGm, SIGh));
ax=gca; ax.XLim(1) = 0;



% obj=findobj('type','line'); for idx=1:length(obj) obj(idx).LineWidth=2; end

