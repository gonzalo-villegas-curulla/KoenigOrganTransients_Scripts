% Extract t10grv, t90grv, t10ft, t90ft

clc, clear;% fprintf('Starting...\n');
load data_proc.mat

% ====== Simulations params ====
pipelist = [1,8,16]; % 1,2,3,...22

% ========= Physical constants =======
rho = 1.2;
co  = 340;          co2 = co^2;
P0  = 820;

% ======= Simulation parameters ==============
fs   = 51.2e3; 
dt   = 1/fs;
Tend = 0.250;
tvec = [0:dt:Tend]';

Nsteps = 10;


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



% Sample select. [1-14,18] (the rest, not valid Spall)(see chart below)
sample_select = 18; 



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

% Effective length corrections *******************************************!!!!!!!!!!!!

Spall_eff = 10^(  -0.03732*12*log2(F1/440) - 4.6202 ); % fitted to computed 
% Spall_eff = data_proc.Spall_eff(sample_select); % Point-wise from measured
% ************************************************************************!!!!!!!!!!!!
Sgrv_eff  = 1*Sgrv_geom; % Couche limite négligeable
Sin_eff   = alpha_in*Sin_geom; % Lambda equation
Sj_eff    = alpha_j*Sj_geom;

% Equation-derived parameters =============================
one_over_A  = data_proc.one_over_Amax(sample_select);
one_over_B  = data_proc.one_over_B(sample_select);
one_over_C  = data_proc.one_over_C(sample_select);
one_over_D  = data_proc.one_over_D(sample_select);
sigMa_full  = data_proc.sigma(sample_select);


results = cell(length(pipelist),1);

% This is the mean values of the given pipe from pipelist[]
tmp = [one_over_A, one_over_B, one_over_C, one_over_D, sigMa_full];
% MXparams = repmat()


%%






% ===================================
% Solve ODE system 
% ===================================
tic

y(1,:) = [1e-5,0];
y(2,:) = [1e-5,0];

tstart = 1e-3;
tfinal = Tend;
refine = 4;

opts  = odeset('RelTol',1e-5,'AbsTol',1e-5,'Refine',refine);


% ===== Solvers:  =====
%   15s, 
%   113 (accurate, sometimes slow)**, 
%   23, 23s
%   avoid 45 (non-adapted time-step)

[t_ode,y] = ode45(@(t_ode,y)...
    solverA(t_ode, y,one_over_A,one_over_B,one_over_C,one_over_D,sigMa_full, ValveRampInit, ValveRampEnd),...
    [tstart tfinal], y(2,:), opts); 
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

pf_over_pgrv_targ = pf(end)/pgrv(end);

% ===================================
%           Plot results
% ===================================
if 0  % Single result plot
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
                        Fact = 0.5 + 0.5*sin(pi*(Tnow-ValveRampInit)/(ValveRampEnd-ValveRampInit) -pi/2);
                    end     
                    ramp(idx) = Fact;
                end
                plot(tax*1e3, ramp, 'linewidth', LW);
                
                legend('Pgrv','Pf','Spall ramp','fontsize',14);
                title(sprintf('PRTgrv= %1.3f ms. PRTf= %1.3f ms',PRTgrv*1e3, PRTf*1e3)   );
                
                % M = [one_over_A,one_over_B;one_over_C,one_over_D];
                % disp(eigs(M))
                % figure; heatmap(M);
end


% ==================================
%            FUNCTIONS 
% ==================================


function dydt = solverA(t_ode, y,one_over_A,one_over_B,one_over_C,one_over_D,sigMa_full, ValveRampInit, ValveRampEnd)

At_factor = one_over_At_val_func(t_ode, ValveRampInit, ValveRampEnd);

one_over_A_time = one_over_A*At_factor;
sigMa_time      = sigMa_full  *At_factor;

dydt = zeros(2,1);
dydt(1) = one_over_A_time*real(sqrt(2*(1-y(1)))) - one_over_B*real(sqrt(y(1)*(2-sigMa_time.^2) - 2*y(2) + sigMa_time.^2)); 
dydt(2) = one_over_C*real(sqrt(y(1)*(2-sigMa_time.^2) - 2*y(2) + sigMa_time.^2)) - one_over_D*real(sqrt(2*y(2)));

end
% -------------------------
function At_fact = one_over_At_val_func(t_ode, ValveRampInit, ValveRampEnd)

    if     t_ode<=ValveRampInit
        At_fact = 0.0;
    elseif ValveRampEnd<t_ode
        At_fact = 1.0;        
    else
        At_fact = 0.5 + 0.5*sin(pi*(t_ode-ValveRampInit)/(ValveRampEnd-ValveRampInit) -pi/2);
    end      

end
