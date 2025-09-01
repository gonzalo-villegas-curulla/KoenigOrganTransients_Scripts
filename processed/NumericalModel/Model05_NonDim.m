%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MODEL 04 (eqs. 9)
%   0. Introduction: 
%      Blah, blah
%
%   1. Constants, variables and parameters:
%   P0 = 820 [Pa] (assumed constant)
%
%   Pgrv, Pf
%   Up2g, Uin, Uj 
% 
%   At , B, C, D, sig 
%
% Gonzalo Villegas Curulla, Tarragona feb 2025 (Running MATLAB R2024b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear;% fprintf('Starting...\n');
load data_proc.mat

% ========= Physical constants =======
rho = 1.2;
co  = 340;          co2 = co^2;
P0  = 820;

% ======= Simulation parameters ==============
fs   = 51.2e2; 
dt   = 1/fs;
Tend = 0.250;
tvec = [0:dt:Tend]';
N    = numel(tvec);

sample_select = 18; % [1-14,18] (the rest, not valid Spall)(see chart below)

% ------------------------------------------------------------------------
%  1   2   3   4   5   6   7   8   9   10   11   12   13   14   15   16     <=== PipeNum
%          1   2   3   4   5       6    7    8         9        10          <=== Sample Num
%          X   X   X   X   X       X    X    X         X         X          <=== Trusted?
% ------------------------------------------------------------------------
%  17   18   19   20   21   22   23   24   25   26   27   28   29   30
%  11        12                       13   14        15        16
%   X         X                        X    X
% ------------------------------------------------------------------------
%  31   32   33   34   35   36   37   38   39   40   41   42   43   44 ...
%       17        18             19        20        21             22
%                  X
% ------------------------------------------------------------------------

% Pallet valve openig time  =============================

ValveRampInit = 0.100; % [s] T-Start opening-ramp pallet valve
ValveRampEnd  = 0.101;  % [s] T-Finish opening-ramp (DROPIC robot time)

% Geometry: retrieve parameters =============================
Spall_geom = data_proc.Spall_geom(sample_select);
Sgrv_geom  = data_proc.Sgrv_geom(sample_select);
Sin_geom   = data_proc.Sin_geom(sample_select);
Sj_geom    = data_proc.Sjet_geom(sample_select);
Vf         = data_proc.Vf(sample_select);
Vgrv       = data_proc.Vgrv(sample_select);
F1         = data_proc.F1(sample_select);


% Geometry: derive parameters =============================

alpha_j  = 1;
alpha_in = alpha_j*(0.5477*Sj_geom/Sin_geom + 0.1398);

% ************************************************************************!!!!!!!!!!!!
        % Spall_eff = 1.0*Spall_geom; %  
Spall_eff = 10^(  -0.03732*12*log2(F1/440) - 4.6202 ); % fitted to computed 
% Spall_eff = data_proc.Spall_eff(sample_select); % Point-wise from measured
% ************************************************************************!!!!!!!!!!!!

Sgrv_eff  = 1*Sgrv_geom;
Sin_eff   = alpha_in*Sin_geom; % Lambda equation
Sj_eff    = alpha_j*Sj_geom;


% Equation-derived parameters =============================

% Ranges
one_over_Amax_range = [5628.05,   14600.6];  % [5628.05,22023.5]
one_over_B_range    = [49.3812,   447.407];  % [19.6542,447.407]
one_over_C_range    = [465.565,   759.252];  % [465.565,1052.01] 
one_over_D_range    = [727.337,   934.094];  % [598.9, 034.094]
Vf_over_Vgrv_range  = [0.0883734, 0.841073]; % in valid region; otherwise [0.0403291,0.841073];
sig_range           = [0.649063,  1.82479];  % in valid region, otherwise [0.649063,2.53989];

one_over_A  = 1.0*Spall_eff*co2/Vgrv*sqrt(rho/P0);            
one_over_B  = 1.0*Sin_eff*co2/Vgrv*sqrt(rho/P0);
one_over_C  = 1.0*Sin_eff*co2/Vf*sqrt(rho/P0);
one_over_D  = 1.0*Sj_eff*co2/Vf*sqrt(rho/P0);
sig_full    = 10.0*Spall_eff/Sgrv_eff; % ******************************?????????????????


% ===================================
% Solve ODE system 
% ===================================
% fprintf('Starting ODE solver...\n');
tic

y(1,:) = [1e-5,0];
y(2,:) = [1e-5,0];

tstart = 1e-3;
tfinal = Tend;
refine = 4;


try 
    opts  = odeset('RelTol',1e-5,'AbsTol',1e-5,'Refine',refine);
    % opts  = odeset('RelTol',1e-5,'AbsTol',1e-5,'Events',@events1,'Events',@events2,'Refine',refine);
    % % % opts  = odeset('RelTol',1e-5,'AbsTol',1e-5,'Events',@events1,'Events',@events2,'Refine',refine,'OutputFcn',@odetpbar);
catch 
    opts  = odeset('RelTol',1e-5,'AbsTol',1e-5,'Events',@events1,'Events',@events2,'Refine',refine);
end



% ===== Solvers:  =====
%   15s, 
%   113 (accurate, sometimes slow)**, 
%   23, 23s
%   avoid 45 (non-adapted time-step)

[t_ode,y] = ode113(@(t_ode,y)...
    solverA(t_ode, y,one_over_A,one_over_B,one_over_C,one_over_D,sig_full, ValveRampInit, ValveRampEnd),...
    [tstart tfinal], y(2,:), opts); 


% fprintf('(done)\n');
toc


yout = y;

% ===================================
%           Analysis
% ===================================

pgrv = yout(:,1);
pf   = yout(:,2);

% Resample homogeneously
pgrv = interp1(t_ode, pgrv, tvec);
pf   = interp1(t_ode, pf, tvec);


t10grv = tvec(find(pgrv/pgrv(end)<0.1,1,'last'));
t90grv = tvec(find(pgrv/pgrv(end)>0.9,1,'first'));
PRTgrv = t90grv-t10grv;

t10f = tvec(find(pf/pf(end)<0.1,1,'last'));
t90f = tvec(find(pf/pf(end)>0.9,1,'first'));
PRTf = t90f-t10f;

% ===================================
%           Plot results
% ===================================

FSZ = 16;
LW  = 1.5;
figure(1);clf;

plot(t_ode*1e3, yout(:,1),'linewidth',LW);
hold on;
plot(t_ode*1e3, yout(:,2),'linewidth',LW);
ylim([-0.125 1.1]);

grid on;
ylabel('Normalized pressure','fontsize', FSZ);
xlabel('Time [ms]', 'FontSize', FSZ);

tax = [0:dt:Tend];
ramp = zeros(length(tax),1);
for idx = 1 : length(tax)
    Tnow = tax(idx);
    if     Tnow<=ValveRampInit
        Fact = 0.0;
    elseif ValveRampEnd<Tnow
        Fact = 1.0;        
    else
%         Fact = 0.5 + 0.5*sin(pi*(Tnow-ValveRampInit)/(ValveRampEnd-ValveRampInit) -pi/2);
        Fact = (Tnow-ValveRampInit)/(ValveRampEnd-ValveRampInit);
%         fprintf('Using a ramp Omega(t) ');
    end     
    ramp(idx) = Fact;
end
plot(tax*1e3, ramp, 'linewidth', LW);

legend('Pgrv','Pf','Spall ramp','fontsize',14);
title(sprintf('PRTgrv= %1.3f ms. PRTf= %1.3f ms',PRTgrv*1e3, PRTf*1e3)   );

M = [one_over_A,one_over_B;one_over_C,one_over_D];
eigs(M)
% figure; heatmap(M);

% ==================================
%            FUNCTIONS 
% ==================================


function dydt = solverA(t_ode, y,one_over_A,one_over_B,one_over_C,one_over_D,sig_full, ValveRampInit, ValveRampEnd)

At_factor = one_over_At_val_func(t_ode, ValveRampInit, ValveRampEnd);

one_over_A_time = At_factor*one_over_A;

%disp(At_factor)
sig_time = sig_full*At_factor;


dydt = zeros(2,1);

dydt(1) = one_over_A_time*real(sqrt(2*(1-y(1)))) - one_over_B*real(sqrt(y(1)*(2-sig_time.^2) - 2*y(2) + sig_time.^2)); 
dydt(2) = one_over_C*real(sqrt(y(1)*(2-sig_time.^2) - 2*y(2) + sig_time.^2)) - one_over_D*real(sqrt(2*y(2)));



end
function At_fact = one_over_At_val_func(t_ode, ValveRampInit, ValveRampEnd)

    if     t_ode<=ValveRampInit
        At_fact = 0.0;
    elseif ValveRampEnd<t_ode
        At_fact = 1.0;        
    else
%         At_fact = 0.5 + 0.5*sin(pi*(t_ode-ValveRampInit)/(ValveRampEnd-ValveRampInit) -pi/2);
        At_fact = (t_ode-ValveRampInit)/(ValveRampEnd-ValveRampInit);
        fprintf('Using a ramp ');
    end      

end

% Zero detection of solution event triggers: % ==========================

function [value,isterminal,direction] = events1(t_ode,y)
value = y(1);     % detect height = 0
isterminal = 1;   % stop the integration
direction = 0;   % negative direction
fprintf('used\n');
end

function [value,isterminal,direction] = events2(t_ode,y)
value = 1e-7;%y(1) | y(2);     % detect height = 0
isterminal = 1;   % stop the integration
direction = 0;   % negative direction
fprintf('used\n');
end
