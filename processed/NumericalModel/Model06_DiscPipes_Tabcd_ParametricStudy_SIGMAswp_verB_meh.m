clc, clear;
load data_proc.mat

% %%%%%%%%%%%%%%%%   CUSTOM USER PARAMETERS  %%%%%%%%%%%%%%%%%%%

pipelist = [8, 18]; % Sampled pipes [1:22]= ["A01",..., "A56"]inventory

% Pallet valve openig time 

ValveRampInit = 0.100;  % [s] T-Start opening-ramp pallet valve
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

NUMVALS = 20;


    sample_select = pipelist(:);

    % PIPE 1 refs

    Ta_ref1 = data_proc.Amax(sample_select(1));
    Tb_ref1 = data_proc.B(sample_select(1));
    Tc_ref1 = data_proc.C(sample_select(1));
    Td_ref1 = data_proc.D(sample_select(1));
    Tomega_ref = 1e-3;

    % Pipe 2 refs
    Ta_ref2 = data_proc.Amax(sample_select(2));
    Tb_ref2 = data_proc.B(sample_select(2));
    Tc_ref2 = data_proc.C(sample_select(2));
    Td_ref2 = data_proc.D(sample_select(2));

    % Parametric slow variation range
    Ta_log     = 1e-3*logspace(-0.15,1.31,NUMVALS);
    Tb_log     = 1e-3*logspace(0.31, 1.7, NUMVALS);
    Tc_lin     = 1e-3*linspace(1.0,2.0, NUMVALS);
    Td_lin     = 1e-3*linspace(1.0, 2.0, NUMVALS);
    Tomega_log = 1e-3*logspace(-1.3, 1.31, 10);

   


    % ===================================
    %           Plot results
    % ===================================

% FUNCTION: [PRTg] = parametrify(Ta, Tb, Tc, Td, Tomega, SIGMA, ValveRampInit, Tend, tvec)    


% ============= Tb ===============


            SIGvec = [0.0 : 0.05 : 1.0]';

fig2 = figure(2); clf; ax2 = axes(fig2); hold on; box on; grid on;
Tb_min     = min(data_proc.B);
Tb_sample8 = data_proc.B(8);

for SIG = [0.0 : 0.05 : 1.0]

    % PIPE 1 (def. 8)
    [PRTg1,~] = parametrifyMULT_func(Ta_ref1, Tb_min, Tc_ref1,Td_ref1, Tomega_ref, SIG, ValveRampInit, Tend, tvec);
    
    try
    plot(ax2, SIG, 1e3*PRTg1, 'ok');
    end    

    
    % PIPE 2 (def. 18)
    
    [PRTg2,~] = parametrifyMULT_func(Ta_ref2, Tb_min, Tc_ref2,Td_ref2, Tomega_ref, SIG, ValveRampInit, Tend, tvec);

    try
    plot(ax2, SIG, 1e3*PRTg2, 'dr');
    end

    drawnow();

end
%



xlabel('$\Sigma$','interpreter','latex');

ylabel('PRT$_g$ [ms]', 'interpreter','latex');
legend({'$A_3^{\#}, \ \Sigma_L$','$A_3^{\#}, \ \Sigma_H$','$A_5^{\#}, \ \Sigma_L $','$A_5^{\#}, \ \Sigma_M $','$A_5^{\#}, \ \Sigma_H$'},...
    'interpreter','latex', 'location','best');
% title(sprintf('SIGMA_{L,M,H}= %1.1f, %1.1f ,%1.1f', SIGl, SIGm, SIGh));
ax=gca; ax.XLim(1) = 0; ax.YLim(1) = 0;



% obj=findobj('type','line'); for idx=1:length(obj) obj(idx).LineWidth=2; end


% =====================================================================================================================
% =====================================================================================================================
% =====================================================================================================================
    


function [PRTg] = parametrify(Ta_pass, Tb_pass, Tc_pass, Td_pass, Tomega_pass, SIG, ValveRampInit, Tend, tvec)

    PRTg = []; 
    PRTg = [PRTg, run_simulation( Ta_pass,Tb_pass,Tc_pass,Td_pass,SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass, tvec)];
        
end


    
function [PRTgrv] = run_simulation(...
       PASS_Amax, ...
       PASS_B, ...
       PASS_C, ...
       PASS_D, ...
       PASS_sigma_full, ...
       Tend, ValveRampInit, ValveRampEnd, tvec)
    
    
    % ===================================
    % Solve ODE system 
    % ===================================
    
    y(1,:) = [1e-5,1.1e-5];
    y(2,:) = [1e-5,1.1e-5];
    
    tstart = 1e-3;
    tfinal = Tend;
    refine = 4;
    
    opts  = odeset('RelTol',1e-5,'AbsTol',1e-5,'Refine',refine);
    opts = odeset(opts,'NonNegative',[1 2]);
    
    % ===== Solvers:  =====
    %   15s, 
    %   113 (accurate, sometimes slow)**, 
    %   23, 23s
    %   45 (non-adapted time-step)
    
    [t_ode,y] = ode45(@(t_ode,y) solverA(t_ode, y,...
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
     
    % Resample homogeneously
    pgrv = interp1(t_ode, pgrv, tvec);
    pf   = interp1(t_ode, pf, tvec);  
    %
    t10grv = tvec(find(pgrv/pgrv(end)<0.1,1,'last'));
    t90grv = tvec(find(pgrv/pgrv(end)>0.9,1,'first'));
    PRTgrv = t90grv-t10grv;
    %
    t10f = tvec(find(pf/pf(end)<0.1,1,'last'));
    t90f = tvec(find(pf/pf(end)>0.9,1,'first'));
    PRTf = t90f-t10f;


    if 0 % Plot all time integrations
                
        figure(10);clf;
        LW = 1.5;
        plot(t_ode*1e3, y(:,1),'linewidth',LW);
        hold on;
        plot(t_ode*1e3, y(:,2),'linewidth',LW);                              
        fprintf(sprintf('Max pgrv: %1.3f. Max pf: %1.3f. Nan grv %d. Nan f %d.\n',...
            max(pgrv), max(pf), isnan(pgrv(end)),isnan(pf(end)) ));
        % xlim([95 120]);
        % ylim([-0.125 1.1]);  
        grid on;
        xlabel('Time [ms]');ylabel('pressure [n.u.]'); 
        title(sprintf('Solution in time. Ta: %1.2f (ms); SIG: %1.1f',PASS_Amax,PASS_sigma_full));
        drawnow();
        % writeVideo(vobj, getframe(gcf) );                
        pause(0.75);
        
    end
   
    
end



% ==================================
%            FUNCTIONS 
% ==================================

% == ODE to solve ======================

function dydt = solverA(t_ode, y,A,B,C,D,sigMa_full, ValveRampInit, ValveRampEnd)

    omeg = omega_func(t_ode, ValveRampInit, ValveRampEnd);    
    dydt = zeros(2,1);
    dydt(1) = omeg*real(sqrt(2*(1-y(1))))/A   -   real(sqrt(y(1)*(2-2*omeg.^2*sigMa_full.^2) -2*y(2) + 2*omeg.^2*sigMa_full^2 ))/B;
    dydt(2) =      real(sqrt(y(1)*(2-2*omeg.^2*sigMa_full.^2) -2*y(2)+2*omeg.^2*sigMa_full.^2 ))/C  -  real(sqrt(2*y(2)))/D;

end

% == Initialisation function Omegat(t) ======================
function OM = omega_func(t_ode, ValveRampInit, ValveRampEnd)

    if     t_ode<=ValveRampInit
        OM = 0.0;
    elseif ValveRampEnd<t_ode
        OM = 1.0;        
    else
        OM = 0.5 - 0.5*cos(pi*(t_ode-ValveRampInit)/(ValveRampEnd-ValveRampInit));
    end      

end