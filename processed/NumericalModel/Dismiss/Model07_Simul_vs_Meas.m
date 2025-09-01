clc, clear;

% %%%%%%%%%%%%%%%%   CUSTOM USER PARAMETERS  %%%%%%%%%%%%%%%%%%%

load data_proc.mat

pipelist = [1:22]; % \in[1,22] (DO NOT INCLUDE MORE THAN 3 PIPES at a time)

% Pallet valve openig time 

ValveRampInit = 0.100; % [s] T-Start opening-ramp pallet valve
ValveRampEnd  = 0.101;  % [s] T-Finish opening-ramp (DROPIC robot time)

% Variation of parameters

lower_bound_factor = 0.8; % Lower end of parameter under modification (Def., 0.5x and 2.0x)
upper_bound_factor = 1.2;

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


% Parameter value allocation and simulation run:

tinit = tic;
results = cell(length(pipelist),1);
MX_results = zeros(length(pipelist), 5); % Ppall, Pgrv, Pf, PRTgrv, PRTf


for pipe_loop_idx = 1:length(pipelist) % ===================================================================================================================


                        % 1/data_proc.Amax(sample_select) ,...
                        % 1/data_proc.B(sample_select),...
                        % 1/data_proc.C(sample_select) ,...
                        % 1/data_proc.D(sample_select) ,...

    sample_select = pipelist(pipe_loop_idx);
    [PRTgrv,PRTf,pf_over_pgrv_targ, Pgrv_trg, Pf_trg, flag_error] = run_simulation( ...
                        data_proc.one_over_Amax(sample_select),...
                        data_proc.one_over_B(sample_select), ...
                        data_proc.one_over_C(sample_select),...
                        data_proc.one_over_D(sample_select),...
                        data_proc.sigma(sample_select),...
                        Tend, ValveRampInit, ValveRampEnd, tvec);
    if flag_error
        MX_results(pipe_loop_idx, :) = nan*ones(1,5);
    else
        MX_results(pipe_loop_idx, :) = [1.0, Pgrv_trg, Pf_trg, PRTgrv, PRTf];        
    end


end 
fprintf('Simulation took %1.3f\n',toc(tinit));


    % ===================================
    %           Plot results
    % ===================================

figure(1); clf;
subplot(2,2,1); % PRTf simul vs meas
errorbar(...
    1e3*flipud(MX_results(:,5)),...
    1e3*data_proc.PRTf_mean,...
    1e3*data_proc.PRTf_std,...
    'v');
grid on; xlabel('Simul [ms]'); ylabel('Meas [ms]'); title('PRT_f');
ylim([0 8]);


subplot(2,2,2); % PRTgrv
errorbar(...
    1e3*flipud(MX_results(:,4)),...
    1e3*data_proc.PRTgrv_mean,...
    1e3*data_proc.PRTgrv_std,...
    'v');
grid on; xlabel('Simul [ms]'); ylabel('Meas [ms]'); title('PRT_{grv}');
ylim([0 8]);


subplot(2,2,3); %Pgrv/Ppall vs simul Pgrv
errorbar(...
    MX_results(:,2),...
    data_proc.Pgrv_mean./data_proc.Ppall_mean,...
    0*data_proc.Pgrv_std./data_proc.Ppall_std,...
    'v');
grid on; xlabel('P_{grv} simul'); ylabel('P_{grv}/P_{pall} meas'); title('P_{grv} ratio');
xlim([0 1.]);ylim([0 1.]);

subplot(2,2,4); % Pf/Pgrv meas vs simul Pf/Pgrv
errorbar(...
    MX_results(:,3)./MX_results(:,2),...
    data_proc.Pf_mean./data_proc.Pgrv_mean,...
    data_proc.Pf_std./data_proc.Pgrv_std,...
    'v');
grid on; xlabel('Simul P_f/P_{grv}'); ylabel('Meas P_f/P_{grv}'); title('P_f/P_{grv} ratio');
hold on;
plot([0,1],[0,1],'-k', 'linewidth',1.5);
xlim([0 1.]);ylim([0 1.3]);
legend('Data','Ideal');

% ====================== ======================

fax = 12*log2(data_proc.F1/440);

figure(2); clf;

% [A]
subplot(2,2,1); % PRTf simul vs meas
errorbar(fax, ...
    1e3*data_proc.PRTf_mean,...
    1e3*data_proc.PRTf_std,...
    'v');
grid on; 
% xlabel('Simul [ms]'); ylabel('Meas [ms]'); 
title('PRT_f');
% ylim([0 8]);
hold on;
plot( fax,...
    1e3*flipud(MX_results(:,5)) );
xlabel('12log_2(F_1/440)');
ylabel('[ms]');
legend('Meas','Simul');

% [B]
subplot(2,2,2); % PRTgrv
errorbar(...
    fax,...
    1e3*data_proc.PRTgrv_mean,...
    1e3*data_proc.PRTgrv_std,...
    'v');
grid on; 
% xlabel('Simul [ms]'); ylabel('Meas [ms]'); 
title('PRT_{grv}');
% ylim([0 8]);
hold on;
plot(fax, 1e3*flipud(MX_results(:,4)) );
xlabel('12log_2(F_1/440)');
ylabel('[ms]');
legend('Meas','Simul');


% [C]
subplot(2,2,3); %Pgrv/Ppall vs simul Pgrv
errorbar(...
    fax,...
    data_proc.Pgrv_mean./data_proc.Ppall_mean,...
    0*data_proc.Pgrv_std./data_proc.Ppall_std,...
    'v');
hold on;
plot(fax, MX_results(:,2));
grid on;  title('P_{grv} ratio');
xlabel('12log_2(F_1/440)');
legend('Meas','Simul');
ylim([0.8 1]);


subplot(2,2,4); % Pf/Pgrv meas vs simul Pf/Pgrv
errorbar(...
    fax,...
    data_proc.Pf_mean./data_proc.Pgrv_mean,...
    data_proc.Pf_std./data_proc.Pgrv_std,...
    'v');
hold on;
plot(fax, MX_results(:,3)./MX_results(:,2))

xlabel('12log_2(F_1/440)'); legend('Meas','Simul');
grid on;
title('P_f/P_{grv} ratio');





% =====================================================================================================================
% =====================================================================================================================
% =====================================================================================================================
    
    
function     [PRTgrv,PRTf,pf_over_pgrv_targ,Pgrv_trg, Pf_trg, flag_error] = run_simulation(...
                       PASS_one_over_Amax,...
                       PASS_one_over_B, ...
                       PASS_one_over_C, ...
                       PASS_one_over_D, ...
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
                                        PASS_one_over_Amax,...
                                        PASS_one_over_B,...
                                        PASS_one_over_C,...
                                        PASS_one_over_D,...
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
    Pgrv_trg = pgrv(end);
    Pf_trg   = pf(end);

    if 0 % Plot on the all all time integrations
                figure(10);clf;
                LW = 1.5;
                plot(t_ode*1e3, yout(:,1),'linewidth',LW);
                hold on;
                plot(t_ode*1e3, yout(:,2),'linewidth',LW);
                ylim([-0.125 1.1]);
                % title(sprintf('f_1 %1.3f', ));
                grid on;
                xlim([0.98*t10grv 1.02*t90f]*1e3);
                drawnow();
                % fprintf(sprintf('Max pgrv: %1.3f. Flag %d. Min pgrv: %1.3f\n',max(pgrv), (max(pgrv)<0), min(pgrv)));
                fprintf('PRTgrv %1.3f [ms], PRtf %1.3f [ms]\n',1e3*PRTgrv, 1e3*PRTf);
                pause();
    end
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
