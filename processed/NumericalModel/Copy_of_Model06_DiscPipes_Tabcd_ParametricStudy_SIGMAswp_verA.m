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

NUMVALS = 10;


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




fig1 = figure(1); clf; ax1 = axes(fig1); hold on; box on; grid on;
Tb_min     = min(data_proc.B);
Tb_sample8 = data_proc.B(8);

for SIG = [0.0 : 0.2 : 1.0]

    % PIPE 1 (def. 8)
    PRTg1a = parametrifyMULT_func(Ta_ref1, Tb_log, Tc_ref1,Td_ref1, Tomega_ref, SIG, ValveRampInit, Tend, tvec);
    
    try
    plot(ax1, 1e3*Tb_log, 1e3*PRTg1a, 'ok');
    end
    xlabel('T_b ms'); ylabel('PRTg ');
    title('PIPE1 & PIPE2: PRTg via parametrifyMULT()');
    drawnow();

    
    % PIPE 2 (def. 18)
    
    [PRTg2a,~] = parametrifyMULT_func(Ta_ref2, Tb_log, Tc_ref2,Td_ref2, Tomega_ref, SIG, ValveRampInit, Tend, tvec);
    PRTg2b     = parametrify(    Ta_ref2, Tb_log, Tc_ref2,Td_ref2, Tomega_ref, SIG, ValveRampInit, Tend, tvec); % HARD CODED BELOW    

    figure();    
    plot(abs(PRTg2a-PRTg2b)./PRTg2a, '-o');grid on;
    title('PIPE2: diff between PRTg via parametrifyMULT() and parametrify()');    
    drawnow();

    try
    plot(ax1, 1e3*Tb_log, 1e3*PRTg2a, 'dr');
    plot(ax1, 1e3*Tb_log, 1e3*PRTg2b, 'dg');
    end    
    drawnow();

end
%

% =====================================================================================================================
% =====================================================================================================================
% =====================================================================================================================
    


function [PRTg] = parametrify(Ta_pass, Tb_pass, Tc_pass, Td_pass, Tomega_pass, SIG, ValveRampInit, Tend, tvec)

    PRTg = []; 
    if length(Ta_pass)>1
        for idx = 1 : length(Ta_pass)
            PRTg = [PRTg, run_simulation( Ta_pass(idx),Tb_pass,Tc_pass,Td_pass,SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass, tvec)];
        end
    end
    if length(Tb_pass)>1
        for idx = 1 : length(Tb_pass)
            PRTg = [PRTg, run_simulation( Ta_pass,Tb_pass(idx),Tc_pass,Td_pass,SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass, tvec)];
        end
    end
    if length(Tc_pass)>1
        for idx = 1 : length(Tc_pass)
            PRTg = [PRTg, run_simulation( Ta_pass,Tb_pass,Tc_pass(idx),Td_pass,SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass, tvec)];
        end
    end
    if length(Td_pass)>1
        for idx = 1 : length(Td_pass)
            PRTg = [PRTg, run_simulation( Ta_pass,Tb_pass,Tc_pass,Td_pass(idx),SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass, tvec)];
        end
    end
    if length(Tomega_pass)>1
        for idx = 1 : length(Tomega_pass)
            PRTg = [PRTg, run_simulation( Ta_pass,Tb_pass,Tc_pass,Td_pass,SIG,     Tend, ValveRampInit, ValveRampInit+Tomega_pass(idx), tvec)];
        end
    end

end

% ===================================================
    
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
    if pgrv(end)<pf(end)
        PRTgrv=NaN;
    end
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
    fprintf('inline func\n');
    omeg = omega_func(t_ode, ValveRampInit, ValveRampEnd);    
    dydt = zeros(2,1);
    dydt(1) = omeg*real(sqrt(2*(1-y(1))))/A   -   real(sqrt(y(1)*(2-2*omeg.^2*sigMa_full.^2) -2*y(2) + 2*omeg.^2*sigMa_full^2 ))/B;
    dydt(2) =      real(sqrt(y(1)*(2-2*omeg.^2*sigMa_full.^2) -2*y(2)+2*omeg.^2*sigMa_full.^2 ))/C  -  real(sqrt(2*y(2)))/D;

end