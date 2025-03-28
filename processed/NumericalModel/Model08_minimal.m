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

% Per form simulation ==================================
for pipe_loop_idx = 1:length(pipelist) 

    sample_select = pipelist(pipe_loop_idx);
    % [PRTgrv,PRTf,pf_over_pgrv_targ, Pgrv_trg, Pf_trg, flag_error] = run_simulation( ...
    %                     1.*data_proc.Amax(sample_select),...
    %                     1.*data_proc.B(sample_select), ...
    %                     data_proc.C(sample_select),...
    %                     data_proc.D(sample_select),...
    %                     data_proc.sigma(sample_select),...
    %                     Tend, ValveRampInit, ValveRampEnd, tvec);

    Xi   = data_proc.Amax(sample_select)/data_proc.B(sample_select);
    Zeta = data_proc.C(sample_select)/data_proc.D(sample_select);
    Volrat = data_proc.B(sample_select)/data_proc.C(sample_select);

    Atmp = 0.8*data_proc.PRTgrv_mean(sample_select); % From PRT
    
    Btmp = Atmp/Xi;    % = data_proc.B(sample_select)./data_proc.Amax(sample_select)*Atmp;
    Ctmp = Btmp/Volrat;% = data_proc.C(sample_select)./data_proc.B(sample_select)*Btmp;
    Dtmp = Ctmp/Zeta;  % = data_proc.D(sample_select)./data_proc.C(sample_select)*Ctmp;

   [PRTgrv,PRTf,pf_over_pgrv_targ, Pgrv_trg, Pf_trg, flag_error] = run_simulation( ...
                        Atmp,...
                        Btmp, ...
                        Ctmp,...
                        Dtmp,...
                        data_proc.sigma(sample_select),...
                        Tend, ValveRampInit, ValveRampEnd, tvec);
   

        MX_results(pipe_loop_idx, :) = [1.0, Pgrv_trg, Pf_trg, PRTgrv, PRTf];   
        res_Pgrv_trg(pipe_loop_idx) = Pgrv_trg;
        res_Pf_trg(pipe_loop_idx)   = Pf_trg;
        res_PRTgrv(pipe_loop_idx)   = PRTgrv;
        res_Pf(pipe_loop_idx)       = PRTf;
        Ares(pipe_loop_idx) = data_proc.Amax(sample_select);
        Bres(pipe_loop_idx) = data_proc.Amax(sample_select);
    
end


    % ===================================
    %           Plot results
    % ===================================

figure(10); clf;
subplot(2,2,1); % PRTf simul vs meas
errorbar(...
    1e3*MX_results(:,5),...
    1e3*data_proc.PRTf_mean,...
    1e3*data_proc.PRTf_std,...
    'v');
grid on; xlabel('Simul [ms]'); ylabel('Meas [ms]'); title('PRT_f');
ylim([0 8]);


subplot(2,2,2); % PRTgrv
errorbar(...
    1e3*MX_results(:,4),...
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

figure(11); clf;

% [A]
subplot(2,2,2); % PRTf simul vs meas
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
    1e3*MX_results(:,5) );
xlabel('12log_2(F_1/440)');
ylabel('[ms]');
legend('Meas','Simul');

% [B]
subplot(2,2,1); % PRTgrv
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
plot(fax, 1e3*MX_results(:,4) );
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
    

% Simulation FUNCTION

function     [PRTgrv,PRTf,pf_over_pgrv_targ,Pgrv_trg, Pf_trg, flag_error] = run_simulation(...
                       PASS_Amax,...
                       PASS_B, ...
                       PASS_C, ...
                       PASS_D, ...
                       PASS_sigma_full, ...
                       Tend, ValveRampInit, ValveRampEnd, tvec)
    
    flag_error = 0;    
    y(1,:) = [1e-5,0];
    y(2,:) = [1e-5,0];
    
    tstart = 1e-3;
    tfinal = Tend;
    refine = 4;
    opts  = odeset('RelTol',1e-5,'AbsTol',1e-5,'Refine',refine);
    
    [t_ode,y] = ode113(@(t_ode,y) solverA(t_ode, y,...
                                        PASS_Amax,...
                                        PASS_B,...
                                        PASS_C,...
                                        PASS_D,...
                                        PASS_sigma_full, ...
                                        ValveRampInit, ValveRampEnd),...
                                        [tstart tfinal], y(2,:), opts); 
    yout = y;
    
            % ====  Analysis  ===============================
            
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
                        figure(20);clf;
                        LW = 1.5;
                        plot(t_ode*1e3, yout(:,1),'linewidth',LW);
                        hold on; grid on;
                        plot(t_ode*1e3, yout(:,2),'linewidth',LW);
                        
                        % xlim([0.98*t10grv 1.02*t90f]*1e3);  
                        xlim([100, 120]);
                        ylim([-0.125 1.1]);
                        title(sprintf('A: %1.2f, B: %1.2f',1e3*PASS_Amax, 1e3*PASS_B));
                        drawnow();
                        pause();
            end
end




% == ODE to solve ======================

function dydt = solverA(t_ode, y, A,B,C,D,sigMa_full, ValveRampInit, ValveRampEnd)

omeg = omega_func(t_ode, ValveRampInit, ValveRampEnd);

dydt = zeros(2,1);
% dydt(1) = omeg*real(sqrt(2*(1-y(1))))/A - real(sqrt(y(1)*(2-omeg.^2*sigMa_full.^2) -2*y(2) + omeg.^2*sigMa_full^2 ))/B;
% dydt(2) = real(sqrt(y(1)*(2-omeg.^2*sigMa_full.^2) -2*y(2)+omeg.^2*sigMa_full.^2 ))/C - real(sqrt(2*y(2)))/D;



% Redo calculation gives a missing factor 2, corrected below:
dydt(1) = omeg*real(sqrt(2*(1-y(1))))/A - real(sqrt(y(1)*(2-2*omeg.^2*sigMa_full.^2) -2*y(2) + 2*omeg.^2*sigMa_full^2 ))/B;
dydt(2) = real(sqrt(y(1)*(2-2*omeg.^2*sigMa_full.^2) -2*y(2)+2*omeg.^2*sigMa_full.^2 ))/C - real(sqrt(2*y(2)))/D;


end
%  == OMEGA function RAMP (0,...,1) ==========================
function OM = omega_func(t_ode, ValveRampInit, ValveRampEnd)

    if     t_ode<=ValveRampInit
        OM = 0.0;
    elseif ValveRampEnd<t_ode
        OM = 1.0;        
    else
        OM = 0.5 + 0.5*sin(pi*(t_ode-ValveRampInit)/(ValveRampEnd-ValveRampInit) -pi/2);
    end      

end
