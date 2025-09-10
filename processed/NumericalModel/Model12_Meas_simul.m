clc, clear;

load data_proc.mat

pipelist = [8]; % \in[1,22]

CHARtime_omega = 7e-3;

% Pallet valve openig time 
ValveRampInit = 0.005; % [s] T-Start opening-ramp pallet valve

RLIMms = 40;

ValveRampEnd_vec = ValveRampInit + CHARtime_omega;


% ========= Physical constants =======
rho = 1.2;
co  = 340;          co2 = co^2;

% ======= Simulation parameters ==============
fs   = 51.2e3; 
dt   = 1/fs;
Tend = 3.300;
tvec = [0:dt:Tend]';


% Simulation ==================================


ValveRampEnd  = ValveRampEnd_vec;  

sample_select = pipelist(1);

Xi   = data_proc.Amax(sample_select)/data_proc.B(sample_select);
Zeta = data_proc.C(sample_select)/data_proc.D(sample_select);
Volrat = data_proc.B(sample_select)/data_proc.C(sample_select);

Atmp = 0.8*data_proc.PRTgrv_mean(sample_select); % From PRT    
Btmp = Atmp/Xi;    % => data_proc.B(sample_select)./data_proc.Amax(sample_select)*Atmp;
Ctmp = Btmp/Volrat;% => data_proc.C(sample_select)./data_proc.B(sample_select)*Btmp;
Dtmp = Ctmp/Zeta;  % => data_proc.D(sample_select)./data_proc.C(sample_select)*Ctmp;

Atmp = data_proc.Amax(sample_select);
Btmp = data_proc.B(sample_select);
Ctmp = data_proc.C(sample_select);
Dtmp = data_proc.D(sample_select);

        PASS_Amax=Atmp;
        PASS_B=Btmp;
        PASS_C=Ctmp;
        PASS_D=Dtmp;
        PASS_sigma_full=data_proc.sigma(sample_select);

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
pgrv_simul = yout(:,1);
pf_simul   = yout(:,2);


% Resample homogeneously
pgrv_simul = interp1(t_ode, pgrv_simul, tvec);
pf_simul   = interp1(t_ode, pf_simul, tvec);  
%
t10grv_simul = tvec(find(pgrv/pgrv(end)<0.1,1,'last'));
t90grv_simul = tvec(find(pgrv/pgrv(end)>0.9,1,'first'));
PRTgrv_simul = t90grv_simul-t10grv_simul;
%
t10f_simul = tvec(find(pf_simul/pf_simul(end) < 0.1,1,'last'));
t90f_simul = tvec(find(pf_simul/pf_simul(end) > 0.9,1,'first'));
PRTf_simul = t90f_simul-t10f_simul;
%
pf_over_pgrv_targ = pf_simul(end)/pgrv_simul(end);
Pgrv_simul_targ         = pgrv_simul(end);
Pf_simul_targ           = pf_simul(end);

omegams = 1e3*(ValveRampEnd_vec-ValveRampInit);

tvec_buffer = tvec;

% =================
% Measured data
% =================
fprintf('Starting to process measured data\n');
currdir = pwd;

PIPENUM = sample_select; TRANSNUM  = 10; 

files = dir('../../A*.mat');
load(fullfile('../../',files(PIPENUM).name), 'PR_data', 'tvec');
tmeas = tvec;

x = PR_data.pressureData;
for idx = 1 : 6
    x(idx,:) = x(idx,:) - mean(x(idx,1:1e3));
end

PRES  = x(2,:); % x(1,:); (Keller)
PGRV  = x(3,:);
PF    = x(4,:);
PRAD  = x(5,:);
KEYV  = x(6,:);

[Lp,Vf,Pw,Tnhd,Sin,Sj,Hm,h,Wm,Rp,palletLHS,palletWid,palletRHS,palletHtraj] = getgeometry(PIPENUM);
Rin     = sqrt(Sin/pi);
Vgroove = Pw*0.51*0.05; % [m^3]
Sslot   = Pw*0.1298;

[KeyDownIdx,KeyUpIdx,KeyMovingTime,DurNotesInS, VelPeakIdxPos,VelPeakIdxNeg] = DetectVelocityPeaks_func(...
    KEYV, tvec, files(PIPENUM).name);


PRECUT      = fix(0.010*fs);
SimulLength = fix(Tend*fs);
ll          = -PRECUT:-PRECUT + SimulLength ; % Mask, selection

pres_meas = PRES(   KeyDownIdx(TRANSNUM)      + ll );
pgrv_meas = PGRV(   KeyDownIdx(TRANSNUM)      + ll );
pf_meas   = PF(     KeyDownIdx(TRANSNUM)      + ll );

% Normalize them by P0 (Ppall)
pgrv_meas = pgrv_meas(:)./pres_meas(:);
pf_meas   = pf_meas(:)  ./pres_meas(:);

% Align them at t10grv and t10f from simul

%"Normalized" target pressures
Pgrv_targ_meas = mean(pgrv_meas(  fix(0.33*length(pgrv_meas))  :  fix(0.66*length(pgrv_meas))   ));
Pf_targ_meas   = mean(pf_meas(   fix(0.33*length(pf_meas))  :  fix(0.66*length(pf_meas))   ));

t10grv_meas_idx = find( pgrv_meas(1:end/2) < 0.1*Pgrv_targ_meas, 1, 'last');
t90grv_meas_idx = find( pgrv_meas(1:end/2) > 0.9*Pgrv_targ_meas, 1, 'first');
PRTgrv_meas     = dt*(t90grv_meas_idx - t10grv_meas_idx);

t10f_meas_idx = find( pf_meas(1:end/2) < 0.1*Pf_targ_meas, 1, 'last');
t90f_meas_idx = find( pf_meas(1:end/2) > 0.9*Pf_targ_meas, 1, 'first');
PRTf_meas     = dt*(t90f_meas_idx-t10f_meas_idx);


pgrv_meas = circshift(pgrv_meas, -t10grv_meas_idx + round(t10grv_simul/dt));
pf_meas   = circshift(pf_meas  , -t10f_meas_idx   + round(t10f_simul/dt));


    % ===================================
    %           Plot results
    % ===================================


figure(1); clf;
plot( (tvec_buffer-t10grv_simul)*1e3, pgrv_meas,'k');  grid on;
title(sprintf('PRTgrv meas %1.3f ms / PRTgrv simul %1.3f ms / OM_{CT}= %1.1f ms',1e3*PRTgrv_meas,1e3*PRTgrv_simul, 1e3*CHARtime_omega));
hold on;
plot((tvec_buffer-t10grv_simul)*1e3, pgrv, '--k');
xlim([-t10grv_simul*1e3 RLIMms]);
xlabel('Time [ms]');
ylabel('Pressure [n.u.]');
legend('P_{grv} measured','P_{grv} simulated', 'location','best');

figure(2); clf;
plot((tvec_buffer-t10f_simul)*1e3, pf_meas,'k'); grid on;
title(sprintf('PRTf meas %1.3f ms/PRTf simul %1.3f ms / OM_{CT} %1.1f ms',1e3*PRTf_meas,1e3*PRTf_simul, 1e3*CHARtime_omega));
hold on;
plot((tvec_buffer-t10f_simul)*1e3, pf, '--k');
xlim([-t10f_simul*1e3 RLIMms]); 
xlabel('Time [ms]');
ylabel('Pressure [n.u]');
legend('P_f measured','P_f simulated','location','best');



















% =====================================================================================================================
% =====================================================================================================================
% =====================================================================================================================
   




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
