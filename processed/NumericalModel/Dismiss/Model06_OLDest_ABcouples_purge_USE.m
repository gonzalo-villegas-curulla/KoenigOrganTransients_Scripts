clc, clear;
load data_proc.mat


% %%%%%%%%%%%%%%%%   CUSTOM USER PARAMETERS  %%%%%%%%%%%%%%%%%%%

pipelist = [8]; % \in[1,22] (DO NOT INCLUDE MORE THAN 3 PIPES at a time)

% Pallet valve openig time 

ValveRampInit = 0.100; % [s] T-Start opening-ramp pallet valve
ValveRampEnd  = 0.101;  % [s] T-Finish opening-ramp (DROPIC robot time)

% Variation of parameters

lower_bound_factor = 0.8; % Lower end of parameter under modification (Def., 0.5x and 2.0x)
upper_bound_factor = 1.2;

SIG = 0.9;

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

rang  = (LOWER : step : UPPER)';

% Parameter value allocation and simulation run:

tinit = tic;
results = cell(length(pipelist),1);
for pipe_loop_idx = 1:length(pipelist) % ===================================================================================================================

    sample_select = pipelist(pipe_loop_idx);
    
    % POPULATE test variables:
    results{pipe_loop_idx}.Amodif.vals     = rang;
    results{pipe_loop_idx}.Bmodif.vals     = rang;   

    % RUN simulations on them:
    for IDX = 1:length(results{pipe_loop_idx}.Amodif.vals)
        for JDX = 1 : length(results{pipe_loop_idx}.Bmodif.vals)
            [PRTgrv,PRTf,pf_over_pgrv_targ, flag_error] = run_simulation( ...
                                results{pipe_loop_idx}.Amodif.vals(IDX) ,...
                                results{pipe_loop_idx}.Bmodif.vals(JDX) ,...
                                data_proc.C(sample_select) ,...
                                data_proc.D(sample_select) ,...
                                SIG,...
                                Tend, ValveRampInit, ValveRampEnd, tvec);
            if flag_error
                MX_PRTgrv(IDX,JDX)  = nan;
                MX_PRTf(IDX,JDX)    = nan;
                MX_Pf_Pgrv(IDX,JDX) = nan;
            else
                MX_PRTgrv(IDX,JDX)  = PRTgrv;
                MX_PRTf(IDX,JDX)    = PRTf;
                MX_Pf_Pgrv(IDX,JDX) = pf_over_pgrv_targ;
            end
        end
        %
    end 
    % Save 
    results{pipe_loop_idx}.results.MX_PRTgrv  = MX_PRTgrv;
    results{pipe_loop_idx}.results.MX_PRTf    = MX_PRTf;
    results{pipe_loop_idx}.results.MX_Pf_Pgrv = MX_Pf_Pgrv;

end 
fprintf('Simulation took %1.3f\n',toc(tinit));


    % ===================================
    %           Plot results
    % ===================================

for idx = 1 : length(results{1}.results.MX_PRTgrv)
    MXtmp = zeros(size(results{1}.results.MX_PRTgrv));
    MXtmp(idx,:) = results{1}.results.MX_PRTgrv(idx,:);
    MXtmp(:,idx) = results{1}.results.MX_PRTgrv(:,idx);
    [aaa,bbb] = max(MXtmp,[],'all','omitnan', 'linear');
    wheremaxgrv(idx) = bbb(end);
end
for idx = 1 : length(results{1}.results.MX_PRTf)
    MXtmp = zeros(size(results{1}.results.MX_PRTf));
    MXtmp(idx,:) = results{1}.results.MX_PRTf(idx,:);
    MXtmp(:,idx) = results{1}.results.MX_PRTf(:,idx);
    [aaa,bbb] = max(MXtmp,[],'all','omitnan','linear');
    wheremaxf(idx) = bbb(end);
end
[AA,BB] = meshgrid(1e3*results{1}.Amodif.vals, 1e3*results{1}.Bmodif.vals);

% SLOPE PRTgrv max / A=B
mask_unique = 1:length(unique(BB(wheremaxgrv)));
slopemaxgrv = BB(wheremaxgrv(mask_unique))./AA(wheremaxgrv(mask_unique));

% SLOPE PRTf max / A=B
mask_unique = 1:length(unique(BB(wheremaxf)));
slopemaxf = BB(wheremaxf(mask_unique))./AA(wheremaxf(mask_unique));



% ============================== B/W
figure(17); clf;

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
plot( min(AAvec,BBvec)+0.07,...
    1e3*MX_PRTf_vec, ...
    '.','color',[1,1,1]*160/255);
legend('$PRT_g$','$PRT_f$','interpreter','latex','autoupdate','off');

datafit1 = polyfit(min(AAvec, BBvec), ...
    1e3*MX_PRTgrv_vec,...
    1);
querypoints = [0.5: 0.5: 20];
fitvals  = polyval(datafit1, querypoints);
hold on;
plot(querypoints, fitvals, '-k');



% ================
MX_PRTf_vec = results{1}.results.MX_PRTf(:);
mask = isnan(MX_PRTf_vec);
MX_PRTf_vec(mask) = [];
AAvec = AA(:); BBvec=BB(:);
AAvec(mask) = []; BBvec(mask) = [];

datafit2 = polyfit( min(AAvec,BBvec),...
    1e3*MX_PRTf_vec, ...
    1);
fitvals2 = polyval(datafit2, querypoints);

plot(querypoints, fitvals2, '-', 'color',[1,1,1]*160/255);

grid on;
xlabel('$min(A,B)\ [ms]$','interpreter','latex');
ylabel('$PRT\ [ms]$','interpreter','latex');
ylim([0 30]);

plot([0,20],[0,20],'--k');





% =====================================================================================================================
% =====================================================================================================================
% =====================================================================================================================
    
    
   function     [PRTgrv,PRTf,pf_over_pgrv_targ, flag_error] = run_simulation(...
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
    
    [t_ode,y] = ode45(@(t_ode,y) solverA(t_ode, y,...
        PASS_Amax,...
        PASS_B, ...
        PASS_C, ...
        PASS_D, ...
        PASS_sigma_full, ...
        ValveRampInit, ValveRampEnd),...
        [tstart tfinal], y(2,:), opts); 
    yout = y;
    
    % ===================================
    %           Analysis
    % ===================================
    
    pgrv = yout(:,1);
    pf   = yout(:,2);
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
    
    pf_over_pgrv_targ = pf(end)/pgrv(end);

    if 1 % Plot on the all all time integrations
                figure(10);clf;
                LW = 1.5;
                plot(t_ode*1e3, yout(:,1),'linewidth',LW);
                hold on;
                plot(t_ode*1e3, yout(:,2),'linewidth',LW);
                ylim([-0.125 1.1]);                
                fprintf(sprintf('Max pgrv: %1.3f. Nan grv %d. Nan f %d. Max pf: %1.3f\n',...
                    max(pgrv), isnan(pgrv(end)),isnan(pf(end)) , max(pf)));
                xlim([95 120]);
                grid on;
                drawnow();
                pause(0.05);
    end
end



% ==================================
%            FUNCTIONS 
% ==================================


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
        OM = 0.5 + 0.5*sin(pi*(t_ode-ValveRampInit)/(ValveRampEnd-ValveRampInit) -pi/2);
    end      

end
