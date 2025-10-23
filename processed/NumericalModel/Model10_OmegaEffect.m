clc, clear;

load data_proc.mat

pipelist = [8]; %Def 8 % \in[1,22]

% Pallet valve openig time 
ValveRampInit = 0.100; % [s] T-Start opening-ramp pallet valve
ValveRampEnd_vec = ValveRampInit +...
    1e-3*[0.05, 0.07,0.1,0.15,0.2,0.25,0.3,...
    0.35, 0.4, 0.45, 0.5,1,2,3,4,5,6,7,8,9,10,...
    12,14,18,22,30,40,50];


% ========= Physical constants =======
rho = 1.2;
co  = 340;          co2 = co^2;

% ======= Simulation parameters ==============
fs   = 4*51.2e3; 
dt   = 1/fs;
Tend = 3.300;
tvec = [0:dt:Tend]';

SIG = 1;


% Parameter value allocation and simulation run:
MX_results = zeros(length(ValveRampEnd_vec), 5); % Ppall, Pgrv, Pf, PRTgrv, PRTf

res_Pgrv_trg = zeros(length(ValveRampEnd_vec), 1);
res_Pf_trg   = zeros(length(ValveRampEnd_vec), 1);
res_PRTgrv   = zeros(length(ValveRampEnd_vec), 1);
res_PRTf     = zeros(length(ValveRampEnd_vec), 1);

% Perform simulation ==================================
for IDX = 1:length(ValveRampEnd_vec) 

    ValveRampEnd  = ValveRampEnd_vec(IDX);  

    sample_select = pipelist(1);
    Pg_hat = data_proc.Pgrv_mean(sample_select)/data_proc.Ppall_mean(sample_select);
    Pf_hat = data_proc.Pf_mean(sample_select)/data_proc.Ppall_mean(sample_select);

    [Ta, Tb, Tc, Td] = generateABCD(SIG, ...
                        Pf_hat, ...
                        Pg_hat, ...
                        data_proc.PRTgrv_mean(sample_select), ...
                        data_proc.Vgrv(       sample_select), ...
                        data_proc.Vf(         sample_select));
    
    

   [PRTgrv,PRTf,pf_over_pgrv_targ, Pgrv_trg, Pf_trg, flag_error] = run_simulation( ...
                        Ta,...
                        Tb, ...
                        Tc,...
                        Td,...
                        SIG,...
                        Tend, ValveRampInit, ValveRampEnd, tvec);
   

        MX_results(IDX,:) = [1.0, Pgrv_trg, Pf_trg, PRTgrv, PRTf];   
        res_Pgrv_trg(IDX) = Pgrv_trg;
        res_Pf_trg(IDX)   = Pf_trg;
        res_PRTgrv(IDX)   = PRTgrv;
        res_PRTf(IDX)     = PRTf;


    
end


    % ===================================
    %           Plot results
    % ===================================

figure(); clf;

omegams = 1e3*(ValveRampEnd_vec-ValveRampInit);

plot(omegams , 1e3*res_PRTgrv);
grid on; hold on; xlabel('$\Omega (t)$ char. rise time [ms]', 'interpreter','latex');
plot(1e3*(ValveRampEnd_vec-ValveRampInit), 1e3*res_PRTf);

ylabel('PRT [ms]');

legend('PRT$_{grv}$','PRT$_f$','location','best', 'Autoupdate','off', 'interpreter','latex'); %axis equal;
ax=gca; 
ax.XLim(1) = 0; 
xlim([0 50]);
ylim([0 25]);
title(sprintf('SIGMA= %1.2f',SIG));



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
                        plot(t_ode*1e3, yout(:,1)/yout(end,1),'linewidth',LW);
                        hold on; grid on;
                        plot(t_ode*1e3, yout(:,2)/yout(end,2),'linewidth',LW);
                        
                        tsegment = tvec(tvec>ValveRampInit & tvec<ValveRampEnd);
                        plot(1e3*tsegment, omega_func(tsegment, ValveRampInit, ValveRampEnd) , 'g');
                        plot( 1e3*tvec(tvec>ValveRampEnd), ones(length(tvec(tvec>ValveRampEnd)),1), 'g');


                        % xlim([0.98*t10grv 1.02*t90f]*1e3);  
                        % xlim([100, 120]);
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
        OM = 0.5 - 0.5*cos(pi*(t_ode-ValveRampInit)/(ValveRampEnd-ValveRampInit));
    end      

end
