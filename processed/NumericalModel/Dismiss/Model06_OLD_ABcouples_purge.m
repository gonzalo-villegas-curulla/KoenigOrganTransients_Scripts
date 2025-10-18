clc, clear;
load data_proc.mat


% %%%%%%%%%%%%%%%%   CUSTOM USER PARAMETERS  %%%%%%%%%%%%%%%%%%%

pipelist = [8];

% Pallet valve openig time 

ValveRampInit = 0.100; % [s] T-Start opening-ramp pallet valve
ValveRampEnd  = 0.101;  % [s] T-Finish opening-ramp (DROPIC robot time)


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

LOWER = 1e-3;
UPPER = 20e-3;
step  = 0.3e-3;

SIG = 1.00;     FIG = 6;

% Parameter value allocation and simulation run:

tinit = tic;
results = cell(length(pipelist),1);
for pipe_loop_idx = 1:length(pipelist) % ===================================================================================================================

    sample_select = pipelist(pipe_loop_idx);

    Pg_hat = data_proc.Pgrv_mean(sample_select)/data_proc.Ppall_mean(sample_select);
    Pf_hat = data_proc.Pf_mean(sample_select)/data_proc.Ppall_mean(sample_select);
    
    [Ta, Tb, Tc, Td] = generateABCD(SIG, ...
                    Pf_hat, ...
                    Pg_hat, ...
                    data_proc.PRTgrv_mean(sample_select), ...
                    data_proc.Vgrv(       sample_select), ...
                    data_proc.Vf(         sample_select));
    
    % POPULATE test variables:


    results{pipe_loop_idx}.Amodif.vals     = [LOWER:step:UPPER]';
    results{pipe_loop_idx}.Bmodif.vals     = [LOWER:step:UPPER]';   



    % RUN simulations on them:
    for IDX = 1:length(results{pipe_loop_idx}.Amodif.vals)
        for JDX = 1 : length(results{pipe_loop_idx}.Bmodif.vals)
            [PRTgrv,PRTf,flag_error] = run_simulation( ...
                                results{pipe_loop_idx}.Amodif.vals(IDX) ,...
                                results{pipe_loop_idx}.Bmodif.vals(JDX) ,...
                                data_proc.C(sample_select) ,...
                                data_proc.D(sample_select) ,...
                                SIG,...
                                Tend, ValveRampInit, ValveRampEnd, tvec);
            if flag_error
                MX_PRTgrv(IDX,JDX)  = nan;
                MX_PRTf(IDX,JDX)    = nan;
            else
                MX_PRTgrv(IDX,JDX)  = PRTgrv;
                MX_PRTf(IDX,JDX)    = PRTf;
            end
        end
        %
    end 
    % Save 
    results{pipe_loop_idx}.results.MX_PRTgrv  = MX_PRTgrv;
    results{pipe_loop_idx}.results.MX_PRTf    = MX_PRTf;

end 
fprintf('Simulation took %1.3f\n',toc(tinit));


    % ===================================
    %           Plot results
    % ===================================

[AA,BB] = meshgrid(1e3*  results{1}.Amodif.vals,...
                   1e3*   results{1}.Bmodif.vals);


% ============================== B/W
figure(FIG); clf;

AAvec = AA(:); BBvec=BB(:);
MX_PRTgrv_vec = results{1}.results.MX_PRTgrv(:);
mask         = isnan(MX_PRTgrv_vec);
MX_PRTgrv_vec(mask) = [];
AAvec(mask) = [];
BBvec(mask) = [];
MX_PRTf_vec = results{1}.results.MX_PRTf(:);
mask = isnan(MX_PRTf_vec);
MX_PRTf_vec(mask) = [];



plot( min(AAvec, BBvec), ...
    1e3*MX_PRTgrv_vec,...
    '.k');
hold on;
plot( max(AAvec, BBvec)+0.07, ...
    1e3*MX_PRTgrv_vec,...
    '.r');
legend({'min(T$_a$,T$_b$)','max(T$_a$,T$_b$)'},'interpreter','latex','autoupdate','off');
% plot( min(AAvec,BBvec)+0.07,...

datafit1 = polyfit(min(AAvec, BBvec), ...
    1e3*MX_PRTgrv_vec,...
    1);
querypoints = [0.5: 0.5: 20];
fitvals  = polyval(datafit1, querypoints);
hold on;
% plot(querypoints, fitvals, '-k');
grid on;
plot([0,20],[0,20],'--k');


xlabel('[ms]','interpreter','latex');
ylabel('PRT$_g$ [ms]','interpreter','latex');
ylim([0 160]);
title(sprintf('SIGMA: %1.3f',SIG));





% =====================================================================================================================
% =====================================================================================================================
% =====================================================================================================================
    
    
   function     [PRTgrv,PRTf, flag_error] = run_simulation(...
       PASS_Amax, ...
       PASS_B, ...
       PASS_C, ...
       PASS_D, ...
       PASS_sigma_full, ...
       Tend, ValveRampInit, ValveRampEnd, tvec)
    
    
    flag_error = 0;
    % ===================================
    % Solve ODE system 
    % ===================================
    
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
    
    [t_ode,y] = ode113(@(t_ode,y) solverA(t_ode, y,...
        PASS_Amax,...
        PASS_B, ...
        PASS_C, ...
        PASS_D, ...
        PASS_sigma_full, ...
        ValveRampInit, ValveRampEnd),...
        [tstart tfinal], y(2,:), opts); 
    
    
    % ===================================
    %           Analysis
    % ===================================
    
    pgrv = y(:,1);
    pf   = y(:,2);
    if max(abs(pf))/max(abs(pgrv))>1  % Parse exp() explosion of solution
        flag_error = 1;
    end
    
    % Resample homogeneously
    pgrv = interp1(t_ode, pgrv, tvec);
    pf   = interp1(t_ode, pf, tvec);  
    
    t10grv = tvec(find(pgrv/pgrv(end)<0.1,1,'last'));
    t90grv = tvec(find(pgrv/pgrv(end)>0.9,1,'first'));
    PRTgrv = t90grv-t10grv;
    
    t10f = tvec(find(pf/pf(end)<0.1,1,'last'));
    t90f = tvec(find(pf/pf(end)>0.9,1,'first'));
    PRTf = t90f-t10f;
   
    
end



% ==================================
%            FUNCTIONS 
% ==================================

% == ODE to solve ======================

function dydt = solverA(t_ode, y,A,B,C,D,sigMa_full, ValveRampInit, ValveRampEnd)

omeg = omega_func(t_ode, ValveRampInit, ValveRampEnd);

dydt = zeros(2,1);
dydt(1) = omeg*real(sqrt(2*(1-y(1))))/A - real(sqrt(y(1)*(2-2*omeg.^2*sigMa_full.^2) -2*y(2) + 2*omeg.^2*sigMa_full^2 ))/B;
dydt(2) = real(sqrt(y(1)*(2-2*omeg.^2*sigMa_full.^2) -2*y(2)+2*omeg.^2*sigMa_full.^2 ))/C - real(sqrt(2*y(2)))/D;


end
% -------------------------
function OM = omega_func(t_ode, ValveRampInit, ValveRampEnd)

    if     t_ode<=ValveRampInit
        OM = 0.0;
    elseif ValveRampEnd<t_ode
        OM = 1.0;        
    else
        OM = 0.5 - 0.5*cos(pi*(t_ode-ValveRampInit)/(ValveRampEnd-ValveRampInit));
    end      

end
