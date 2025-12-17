clc, clear;
load data_proc.mat

pipelist      = [1:22];
ValveRampInit = 0.100; % [s]
ValveRampEnd  = 0.101;  % [s]


% ========= Physical constants =======
rho = 1.2;
co  = 340;          co2 = co^2;

% ======= Simulation parameters ==============
fs   = 51.2e3; 
dt   = 1/fs;
Tend = 0.500;
tvec = [0:dt:Tend];

simul_all_PRTgrv   = zeros(length(pipelist), 1);
simul_all_PRTf     = zeros(length(pipelist), 1);

SIG = 1.0;

% Perform simulation ==================================
for pipe_loop_idx = 1:length(pipelist) 

    sample_select = pipelist(pipe_loop_idx);

    Pg_hat = data_proc.Pgrv_mean(sample_select)/data_proc.Ppall_mean(sample_select);
    Pf_hat = data_proc.Pf_mean(sample_select)/data_proc.Ppall_mean(sample_select);
    
    [Ta, Tb, Tc, Td] = generateABCD(SIG, ...
                        Pf_hat, ...
                        Pg_hat, ...
                        data_proc.PRTgrv_mean(sample_select), ...
                        data_proc.Vgrv(       sample_select), ...
                        data_proc.Vf(         sample_select));


   [PRTgrv,PRTf] = run_simulation( ...
                        Ta,...
                        Tb, ...
                        Tc,...
                        Td,...
                        SIG,...
                        Tend, ValveRampInit, ValveRampEnd, tvec);
   
    simul_all_PRTgrv(  pipe_loop_idx)   = PRTgrv;
    simul_all_PRTf(    pipe_loop_idx)   = PRTf;
        
    
end


    % ===================================
    %           Plot results
    % ===================================

fax = 12*log2(data_proc.F1/440);

figure(1); clf;

plot(fax, ...
    1e3*data_proc.PRTgrv_mean,...
    '*k');
grid on; hold on;
plot( fax,...
    1e3*simul_all_PRTgrv, ...
    'dk', 'markerfacecolor','k');
xlabel('$12\times log_2(F_1/440)$', 'interpreter','latex');
ylabel('[ms]', 'interpreter','latex');
ylim([0 13]);
legend('PRT$_g$ Meas.','PRT$_g$ Simul.', 'interpreter','latex');

%
figure(2); clf;
plot(...
    fax,...
    1e3*data_proc.PRTf_mean,...    
    '*k');
grid on; hold on;
plot(fax, 1e3*simul_all_PRTf,...
    'dk', 'markerfacecolor','k');
xlabel('$12\times log_2(F_1/440)$', 'interpreter','latex');
ylabel('[ms]', 'interpreter','latex');
ylim([0 13]);
legend('PRT$_f$ Meas.','PRT$_f$ Simul.','interpreter','latex');






% =====================================================================================================================
% =====================================================================================================================
% =====================================================================================================================
    

% Simulation FUNCTION

function     [PRTgrv,PRTf] = run_simulation(...
                       PASS_Amax,...
                       PASS_B, ...
                       PASS_C, ...
                       PASS_D, ...
                       PASS_sigma_full, ...
                       Tend, ValveRampInit, ValveRampEnd, tvec)
     
    y(1,:) = [1e-5,0]; % Initialisae 1st two time steps
    y(2,:) = [1e-5,0]; % Init 1st two time steps
    
    tstart = 1e-6; 
    tfinal = Tend;
    refine = 4;
    opts   = odeset('RelTol',1e-8,'AbsTol',1e-8,'Refine',refine);
    
    [t_ode,y] = ode78(@(t_ode,y) solverA(t_ode, y,...
                                        PASS_Amax,...
                                        PASS_B,...
                                        PASS_C,...
                                        PASS_D,...
                                        PASS_sigma_full, ...
                                        ValveRampInit, ValveRampEnd),...
                                        [tstart tfinal], y(2,:), opts); 
    
    
            % ====  Analysis  ===============================
            
            pgrv = y(:,1);
            pf   = y(:,2);
            
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

function dydt = solverA(t_ode, y, A,B,C,D,sigMa_full, ValveRampInit, ValveRampEnd)

    omeg = omega_func(t_ode, ValveRampInit, ValveRampEnd);
    
    % y(1): Pgrv
    % y(2): Pf
    dydt = zeros(2,1);    
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
        OM = 0.5 - 0.5*cos(pi*(t_ode-ValveRampInit)/(ValveRampEnd-ValveRampInit));
    end      

end
