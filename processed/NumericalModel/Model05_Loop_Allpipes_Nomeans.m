clc, clear;
load data_proc.mat


% %%%%%%%%%%%%%%%%   CUSTOM USER PARAMETERS  %%%%%%%%%%%%%%%%%%%

pipelist = [2,7,14]; % \in[1,22] (DO NOT INCLUDE MORE THAN 3 PIPES at a time)

% Pallet valve openig time 

ValveRampInit = 0.100; % [s] T-Start opening-ramp pallet valve
ValveRampEnd  = 0.101;  % [s] T-Finish opening-ramp (DROPIC robot time)

% Variation of parameters

lower_bound_factor = 0.01; % Lower end of parameter under modification (Def., 0.5x and 2.0x)
upper_bound_factor = 16.0;

Nsteps = 100; % Between min and max range of the parameter variation

flag_means = true;

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

results = cell(length(pipelist),1);
for pipe_loop_idx = 1:length(pipelist) % ===================================================================================================================

    sample_select = pipelist(pipe_loop_idx);
    
    results{pipe_loop_idx}.Amodif.vals = linspace(data_proc.one_over_Amax(sample_select)*lower_bound_factor,...
        data_proc.one_over_Amax(sample_select)*upper_bound_factor,...
        Nsteps);

    results{pipe_loop_idx}.Bmodif.vals = linspace(data_proc.one_over_B(sample_select)*lower_bound_factor,...
        data_proc.one_over_B(sample_select)*upper_bound_factor,...
        Nsteps);

    results{pipe_loop_idx}.Cmodif.vals = linspace(data_proc.one_over_C(sample_select)*lower_bound_factor,...
        data_proc.one_over_C(sample_select)*upper_bound_factor,...
        Nsteps);

    results{pipe_loop_idx}.Dmodif.vals = linspace(data_proc.one_over_D(sample_select)*lower_bound_factor,...
        data_proc.one_over_D(sample_select)*upper_bound_factor,...
        Nsteps);

    results{pipe_loop_idx}.Sigmamodif.vals = linspace(data_proc.sigma(sample_select)*lower_bound_factor,...
        data_proc.sigma(sample_select)*upper_bound_factor,...
        Nsteps);
   

    for IDX = 1:Nsteps
        % Run variation of A on pipe pipelist(pipe_loop_idx)
        [PRTgrv,PRTf,pf_over_pgrv_targ, flag_error] = run_simulation( results{pipe_loop_idx}.Amodif.vals(IDX) ,...
                            data_proc.one_over_B(sample_select),data_proc.one_over_C(sample_select),data_proc.one_over_D(sample_select),data_proc.sigma(sample_select), Tend, ValveRampInit, ValveRampEnd, tvec);
        if flag_error
            MX_PRTgrv(IDX,1)  = nan;
            MX_PRTf(IDX,1)    = nan;
            MX_Pf_Pgrv(IDX,1) = nan;
        else
            MX_PRTgrv(IDX,1)  = PRTgrv;
            MX_PRTf(IDX,1)    = PRTf;
            MX_Pf_Pgrv(IDX,1) = pf_over_pgrv_targ;
        end
        % Run variation of B        
        [PRTgrv,PRTf,pf_over_pgrv_targ, flag_error] = run_simulation(data_proc.one_over_Amax(sample_select),...
            results{pipe_loop_idx}.Bmodif.vals(IDX),...
            data_proc.one_over_C(sample_select),data_proc.one_over_D(sample_select),data_proc.sigma(sample_select) , Tend, ValveRampInit, ValveRampEnd, tvec);
        if flag_error
            MX_PRTgrv(IDX,2)  = nan;
            MX_PRTf(IDX,2)    = nan;
            MX_Pf_Pgrv(IDX,2) = nan;
        else
            MX_PRTgrv(IDX,2)  = PRTgrv;
            MX_PRTf(IDX,2)    = PRTf;
            MX_Pf_Pgrv(IDX,2) = pf_over_pgrv_targ;
        end
        % Run variation of C
        [PRTgrv,PRTf,pf_over_pgrv_targ, flag_error] = run_simulation(data_proc.one_over_Amax(sample_select), data_proc.one_over_B(sample_select),...
                results{pipe_loop_idx}.Cmodif.vals(IDX),...
                data_proc.one_over_D(sample_select), data_proc.sigma(sample_select) , Tend, ValveRampInit, ValveRampEnd, tvec);
        if flag_error
            MX_PRTgrv(IDX,3)  = nan;
            MX_PRTf(IDX,3)    = nan;
            MX_Pf_Pgrv(IDX,3) = nan;
        else
            MX_PRTgrv(IDX,3)  = PRTgrv;
            MX_PRTf(IDX,3)    = PRTf;
            MX_Pf_Pgrv(IDX,3) = pf_over_pgrv_targ;
        end
        % Run variation of D
        [PRTgrv,PRTf,pf_over_pgrv_targ, flag_error] = run_simulation(data_proc.one_over_Amax(sample_select), data_proc.one_over_B(sample_select),data_proc.one_over_C(sample_select), ...        
                results{pipe_loop_idx}.Dmodif.vals(IDX),...                
                data_proc.sigma(sample_select) , Tend, ValveRampInit, ValveRampEnd, tvec);
        if flag_error
            MX_PRTgrv(IDX,4)  = nan;
            MX_PRTf(IDX,4)    = nan;
            MX_Pf_Pgrv(IDX,4) = nan;
        else
            MX_PRTgrv(IDX,4)  = PRTgrv;
            MX_PRTf(IDX,4)    = PRTf;
            MX_Pf_Pgrv(IDX,4) = pf_over_pgrv_targ;
        end
        % Run variation of Sigma
        [PRTgrv,PRTf,pf_over_pgrv_targ, flag_error] = run_simulation(data_proc.one_over_Amax(sample_select), data_proc.one_over_B(sample_select),data_proc.one_over_C(sample_select), data_proc.one_over_D(sample_select),...                
                results{pipe_loop_idx}.Sigmamodif.vals(IDX) , Tend, ValveRampInit, ValveRampEnd, tvec);
        if flag_error
            MX_PRTgrv(IDX,5)  = nan;
            MX_PRTf(IDX,5)    = nan;
            MX_Pf_Pgrv(IDX,5) = nan;
        else
            MX_PRTgrv(IDX,5)  = PRTgrv;
            MX_PRTf(IDX,5)    = PRTf;
            MX_Pf_Pgrv(IDX,5) = pf_over_pgrv_targ;
        end
        %
    end 
    % Save ABCDsigma variational results for pipe pipelist(pipe_loop_idx)
    results{pipe_loop_idx}.results.MX_PRTgrv  = MX_PRTgrv;
    results{pipe_loop_idx}.results.MX_PRTf    = MX_PRTf;
    results{pipe_loop_idx}.results.MX_Pf_Pgrv = MX_Pf_Pgrv;

end 


    % ===================================
    %           Plot results
    % ===================================

% ===============================================================

lgd = cell(length(pipelist), 1);
for idx = 1 : length(pipelist)
    lgd{idx} = sprintf('Pipe %d', pipelist(idx));
end

% PARAMETER =====  A ========
figure(1); clf; 
%
axhA(1) = subplot(1,3,1); hold on;
for idx=1:length(pipelist)
    plot( 1./(results{idx}.Amodif.vals),  1e3*results{idx}.results.MX_PRTgrv(:,1)  );
end
if flag_means
    for idx = 1:length(pipelist)
        plot( 1./data_proc.one_over_Amax(pipelist(idx)) ,...
        1e3*interp1(1./(results{idx}.Amodif.vals),results{idx}.results.MX_PRTgrv(:,1),1./data_proc.one_over_Amax(pipelist(idx)), 'pchip'),...
        'dk');
    end
end
xlabel('A [s]'); ylabel('[ms]'); title('PRTgrv'); grid on; box on;
legend(lgd);
% -----------------------------
        axhA(2) = subplot(1,3,2); hold on;
        for idx=1:length(pipelist)
            plot( 1./(results{idx}.Amodif.vals),  1e3*results{idx}.results.MX_PRTf(:,1)  );
        end        
        if flag_means
             for idx = 1:length(pipelist)
                plot( 1./data_proc.one_over_Amax(pipelist(idx)) ,...
                1e3*interp1(1./(results{idx}.Amodif.vals),results{idx}.results.MX_PRTf(:,1),1./data_proc.one_over_Amax(pipelist(idx)), 'pchip'),...
                'dk');
             end
        end
        xlabel('A [s]'); ylabel('[ms]'); title('PRTf');grid on; box on;
        legend(lgd);
% -----------------------------
axhA(3) = subplot(1,3,3); hold on;
for idx=1:length(pipelist)
    plot( 1./(results{idx}.Amodif.vals),  results{idx}.results.MX_Pf_Pgrv(:,1)  );
end
if flag_means
    for idx = 1:length(pipelist)
        plot( 1./data_proc.one_over_Amax(pipelist(idx)) ,...
        1e3*interp1(1./(results{idx}.Amodif.vals),results{idx}.results.MX_Pf_Pgrv(:,1),1./data_proc.one_over_Amax(pipelist(idx)), 'pchip'),...
        'dk');
    end
end

xlabel('A [s]'); title('Pf/Pgrvf Targ'); ylim([0 1]); grid on; box on;
legend(lgd);
linkaxes(axhA,'x');

% PARAMETER =====  B ========
figure(2); clf; 
%
axhB(1) = subplot(1,3,1); hold on;
for idx=1:length(pipelist)
    plot( 1./(results{idx}.Bmodif.vals),  1e3*results{idx}.results.MX_PRTgrv(:,2)  );
end
xlabel('B [s]'); ylabel('[ms]'); title('PRTgrv'); ylim([0 5.5]);
grid on; box on;
legend(lgd);
% -----------------------------
        axhB(2) = subplot(1,3,2); hold on;
        for idx=1:length(pipelist)
            plot( 1./(results{idx}.Bmodif.vals),  1e3*results{idx}.results.MX_PRTf(:,2)  );
        end
        xlabel('B [s]'); ylabel('[ms]'); title('PRTf'); grid on; box on; ylim([0 5.5]);
        legend(lgd);
% -----------------------------
axhB(3) = subplot(1,3,3); hold on;
for idx=1:length(pipelist)
    plot( 1./(results{idx}.Bmodif.vals),  results{idx}.results.MX_Pf_Pgrv(:,2)  );
end
xlabel('B [s]'); title('Pf/Pgrv Targ'); grid on; box on; ylim([0 1]);
legend(lgd);
linkaxes(axhB,'x');




% PARAMETER =====  C ========
figure(3); clf; 
%
axhC(1) = subplot(1,3,1); hold on;
for idx=1:length(pipelist)
    plot( 1./(results{idx}.Cmodif.vals),  1e3*results{idx}.results.MX_PRTgrv(:,3)  );
end
xlabel('C [s]'); ylabel('[ms]'); title('PRTgrv'); grid on; box on; ylim([0 5.5]);
legend(lgd);
% -----------------------------
        axhC(2) = subplot(1,3,2); hold on;
        for idx=1:length(pipelist)
            plot( 1./(results{idx}.Cmodif.vals),  1e3*results{idx}.results.MX_PRTf(:,3)  );
        end
        xlabel('C [s]'); ylabel('[ms]'); title('PRTf'); grid on; box on; ylim([0 5.5]);
        legend(lgd);
% -----------------------------
axhC(3) = subplot(1,3,3); hold on;
for idx=1:length(pipelist)
    plot( 1./(results{idx}.Cmodif.vals),  results{idx}.results.MX_Pf_Pgrv(:,3)  );
end
xlabel('C[s] (log_{10})'); title('Pf/Pgrv Targ'); ylabel('(log_{10})'); grid on; box on;
legend(lgd);
linkaxes(axhC([1,2]),'x');


% PARAMETER =====  D ========
figure(4); clf; 
%
axhD(1) = subplot(1,3,1); hold on;
for idx=1:length(pipelist)
    plot( 1./(results{idx}.Dmodif.vals),  1e3*results{idx}.results.MX_PRTgrv(:,4)  );
end
xlabel('D [s]'); ylabel('[ms]'); title('PRTgrv'); grid on; box on; ylim([0 5.5]);
legend(lgd);
% -----------------------------
        axhD(2) = subplot(1,3,2); hold on;
        for idx=1:length(pipelist)
            plot( 1./(results{idx}.Dmodif.vals),  1e3*results{idx}.results.MX_PRTf(:,4)  );
        end
        xlabel('D[s]'); ylabel('[ms]'); title('PRTf'); grid on; box on; ylim([0 5.5]);
        legend(lgd);
% -----------------------------
axhD(3) = subplot(1,3,3); hold on;
for idx=1:length(pipelist)
    plot( 1./(results{idx}.Dmodif.vals),  results{idx}.results.MX_Pf_Pgrv(:,4)  );
end
xlabel('D [s] (log_{10})'); ylabel('(log_{10})'); title('Pf/Pgrv Targ'); grid on; box on;
legend(lgd);
linkaxes(axhD([1,2]),'x');


% PARAMETER =====  Sigma ========
figure(5); clf; 
%
subplot(1,3,1); hold on;
for idx=1:length(pipelist)
    plot( (results{idx}.Sigmamodif.vals),  1e3*results{idx}.results.MX_PRTgrv(:,5)  );
end
xlabel('\Sigma'); ylabel('[ms]'); title('PRTgrv'); grid on; box on; ylim([0 5.5]);
legend(lgd);
% -----------------------------
        subplot(1,3,2); hold on;
        for idx=1:length(pipelist)
            plot( results{idx}.Sigmamodif.vals,  1e3*results{idx}.results.MX_PRTf(:,5)  );
        end
        xlabel('\Sigma'); ylabel('[ms]'); title('PRTf'); grid on; box on; ylim([0 5.5]);
        legend(lgd);
% -----------------------------
subplot(1,3,3); hold on;
for idx=1:length(pipelist)
    plot( results{idx}.Sigmamodif.vals,  results{idx}.results.MX_Pf_Pgrv(:,5)  );
end
xlabel('\Sigma'); title('Pf/Pgrv Targ'); grid on; box on; ylim([0 1]);
legend(lgd);



% =====================================================================================================================
% =====================================================================================================================
% =====================================================================================================================
    
    
   function     [PRTgrv,PRTf,pf_over_pgrv_targ, flag_error] = run_simulation(PASS_one_over_Amax, PASS_one_over_B, PASS_one_over_C, PASS_one_over_D, PASS_sigma_full, Tend, ValveRampInit, ValveRampEnd, tvec)
    
    
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
    
    [t_ode,y] = ode45(@(t_ode,y)...
        solverA(t_ode, y,PASS_one_over_Amax,PASS_one_over_B,PASS_one_over_C,PASS_one_over_D,PASS_sigma_full, ValveRampInit, ValveRampEnd),...
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

    if 0 % Plot on the all all time integrations
                figure(1);clf;
                LW = 1.5;
                plot(t_ode*1e3, yout(:,1),'linewidth',LW);
                hold on;
                plot(t_ode*1e3, yout(:,2),'linewidth',LW);
                ylim([-0.125 1.1]);
                drawnow();
                fprintf(sprintf('Max pgrv: %1.3f. Flag %d. Min pgrv: %1.3f\n',max(pgrv), (max(pgrv)<0), min(pgrv)));
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
