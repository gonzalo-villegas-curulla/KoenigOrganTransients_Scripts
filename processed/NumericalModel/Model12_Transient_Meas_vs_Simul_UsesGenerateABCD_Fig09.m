clc;
clearvars -except MX1 sample_select TRANSNUM IIDX JDX
load data_proc.mat;
% % % % % % sample_select = 8; % A#_0 % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
CHARtime_omega = 1e-3;

% Pallet valve openig time 
ValveRampInit = 0.010; % [s] T-Start opening-ramp pallet valve
RLIMms        = 20;

ValveRampEnd  = ValveRampInit + CHARtime_omega;


% ========= Physical constants =======
rho = 1.2;
co  = 340;          co2 = co^2;

% ======= Simulation parameters ==============
fs   = 51.2e3; 
dt   = 1/fs;
Tend = 3.300; 
% tmeas = [0:dt:Tend]';

SIGvec = [0.0, 1.0]';
SIGvec = [1.0]';% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%<<<<<<<<<<<<<<<<<<<<<


% ========== Figure handles =============
fh1 = figure(1); clf; axh1 = axes(fh1); hold on; grid on; box on;
fh2 = figure(2); clf; axh2 = axes(fh2); hold on; grid on; box on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%               %%%%%%%%%%%%%%%%%%%%%%%%%%
%               %     Measured data      %
%               %%%%%%%%%%%%%%%%%%%%%%%%%%


Pg_hat = data_proc.Pgrv_mean(sample_select)/data_proc.Ppall_mean(sample_select);
Pf_hat = data_proc.Pf_mean(sample_select)/data_proc.Ppall_mean(sample_select);

% =================
% Measured data
% =================
fprintf('Starting to process measured data\n');



PIPENUM = sample_select; 
% % % % % % TRANSNUM  = 10; %Def: 10  %%%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<





files   = dir('../../A*.mat');
load(fullfile('../../',files(PIPENUM).name), 'PR_data', 'tvec');

x = PR_data.pressureData;
for idx = 1 : 6
    x(idx,:) = x(idx,:) - mean(x(idx,1:1e3));
end


% x(1,:); (Keller)
PRES  = x(2,:); % Endevco1
PGRV  = x(3,:); % Endevco2
PF    = x(4,:); % Endevco3
PRAD  = x(5,:); % B&K
KEYV  = x(6,:); % Polytec Vibr PVD100

[Lp,Vf,Pw,Tnhd,Sin,Sj,Hm,h,Wm,Rp,palletLHS,palletWid,palletRHS,palletHtraj] = getgeometry(PIPENUM);
Rin     = sqrt(Sin/pi);
Vgroove = Pw*0.51*0.05; % [m^3]
Sslot   = Pw*0.1298;

[KeyDownIdx,KeyUpIdx,KeyMovingTime,DurNotesInS, VelPeakIdxPos,VelPeakIdxNeg] = DetectVelocityPeaks_func(...
    KEYV, tvec, files(PIPENUM).name);


PRECUT      = fix(0.010*fs);
SimulLength = fix(Tend*fs);
ll          = -PRECUT:-PRECUT + SimulLength ; 

pres_meas = PRES(   KeyDownIdx(TRANSNUM)      + ll );
pgrv_meas = PGRV(   KeyDownIdx(TRANSNUM)      + ll );
pf_meas   = PF(     KeyDownIdx(TRANSNUM)      + ll );

L = length(pgrv_meas);
tmeas = dt*[0:L-1];

% Normalize them by P0 (Ppall)
pgrv_meas = pgrv_meas(:)./pres_meas(:);
pf_meas   = pf_meas(:)  ./pres_meas(:);

pgrv_meas_targ = mean(pgrv_meas(fix(0.33*length(pgrv_meas)):fix(0.66*length(pgrv_meas))));
tmeas_p = tmeas - dt*find(pgrv_meas(1:fix(0.5*length(pgrv_meas)))/pgrv_meas_targ<0.1,1, 'last');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%               %%%%%%%%%%%%%%%%%%%%%%%%%%
%               %    Simulated data      %
%               %%%%%%%%%%%%%%%%%%%%%%%%%%


plot(axh1, tmeas_p*1e3, pgrv_meas, 'k');
plot(axh2, tmeas_p*1e3, pf_meas, 'k');


for IDX = 1 : length(SIGvec) 

    SIG = SIGvec(IDX);
    [Ta, Tb, Tc, Td] = generateABCD(SIG, ...
        Pf_hat, ...
        Pg_hat, ...
        data_proc.PRTgrv_mean(sample_select), ...
        data_proc.Vgrv(sample_select), ...
        data_proc.Vf(sample_select) );
    
    % Simu: initial conditions
    y(1,:) = [1e-5,0];
    y(2,:) = [1e-5,0];
    tstart = 1e-3;
    tfinal = Tend;
    refine = 4;
    opts   = odeset('RelTol',1e-5,'AbsTol',1e-5,'Refine',refine);
    
    [tsimu,y] = ode113(@(t_ode,y) solverA(t_ode, y,...
                                    Ta,...
                                    Tb,...
                                    Tc,...
                                    Td,...
                                    SIG, ...
                                    ValveRampInit, ValveRampEnd),...
                                    [tstart tfinal], y(2,:), opts); 

    % ====  Analysis  ===============================
    pgrv_simu = y(:,1);
    pf_simu   = y(:,2);
    
    
    % Homogeneously sample data 
    pgrv_simu = interp1(tsimu, pgrv_simu, tmeas);
    pf_simu   = interp1(tsimu, pf_simu,  tmeas);    
    tsimu     = interp1(tsimu, tsimu, tmeas);
    
    
    t10grv_simu = tsimu(find(pgrv_simu/max(pgrv_simu) < 0.1,1,'last'));
    tsimu_p     = tsimu - 0*t10grv_simu;
    
    
        % ===================================
        %           Plot simul
    everyN = 100; % Marker density
    
    
    
    plot(axh1, tsimu_p*1e3, pgrv_simu, ...
        '--', 'Marker','o','MarkerIndices',1:everyN:length(pgrv_simu));
    
    plot(axh2, tsimu_p *1e3, pf_simu,...
          '--', 'Marker','o',...
          'MarkerIndices',1:everyN:length(pf_simu) );

    % ====================
    % Sigma = 0,1

    
    finit = zeros(size(tsimu));
    for idx = 1 : length(finit)
        finit(idx) = omega_func(tsimu(idx), ValveRampInit, ValveRampEnd);
    end
    finit(isnan(finit)) = 0;
    finit(end) = 1;


    t10key_simu = tsimu(find(finit<0.1, 1, 'last')); 
    t10g_simu   = t10grv_simu; 
    t10f_simu   = tsimu(find( pf_simu/max(pf_simu)<0.1, 1, 'last' ));
    % taumech_simu = t10g_simu - t10key_simu;
    tauflow_simu = t10f_simu - t10g_simu;
    


end

XLIM = [-5 10];
XLIM = [-5 30];
xlim(axh1, XLIM);
xlim(axh2, XLIM);


ylim(axh1, [-0.05  1.0]);
ylim(axh2, [-0.05, 0.4]);


xlabel(axh1, 'Time [ms]', 'interpreter','latex');
ylabel(axh1, '$\overline{P}_g \ [n.u.]$', 'interpreter','latex');
legend(axh1, {'Measured','Sim. $\Sigma = 0$','Sim. $\Sigma = 1$'}, ...
    'location','southeast','interpreter','latex');




xlabel(axh2, 'Time [ms]', 'interpreter','latex');
ylabel(axh2, '$\overline{P}_f \ [n.u]$', 'interpreter','latex');
legend(axh2, {'Measured','Sim. $\Sigma = 0$','Sim. $\Sigma = 1$'},...
    'location','southeast','interpreter','latex');
drawnow();

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
