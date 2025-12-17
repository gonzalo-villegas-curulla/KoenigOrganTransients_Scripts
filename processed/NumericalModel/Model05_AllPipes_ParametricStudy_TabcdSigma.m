clc, clear;
load data_proc.mat


% %%%%%%%%%%%%%%%%   CUSTOM USER PARAMETERS  %%%%%%%%%%%%%%%%%%%

pipelist = [1,8,12]; % \in[1,22] (DO NOT INCLUDE MORE THAN 3 PIPES at a time)

YLIMS = [0 5.5];

% Pallet valve openig time 

ValveRampInit = 0.100;  % [s] T-Start opening-ramp pallet valve
ValveRampEnd  = 0.101;  % [s] T-Finish opening-ramp (DROPIC robot time)

% Variation of parameters

lower_bound_factor = 0.8; % Lower end of parameter under modification (Def., 0.5x and 2.0x)
upper_bound_factor = 1.2;

Nsteps = 10; % Between min and max range of the parameter variation

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


tinit = tic;
results = cell(length(pipelist),1);
for pipe_loop_idx = 1:length(pipelist) % ===================================================================================================================

    sample_select = pipelist(pipe_loop_idx);
    
    results{pipe_loop_idx}.Amodif.vals     = linspace(min(data_proc.Amax)*lower_bound_factor, max(data_proc.Amax)*upper_bound_factor, Nsteps)';
    results{pipe_loop_idx}.Bmodif.vals     = linspace(min(data_proc.B)*lower_bound_factor, max(data_proc.B)*upper_bound_factor, Nsteps)';
    results{pipe_loop_idx}.Cmodif.vals     = linspace(min(data_proc.C)*lower_bound_factor, max(data_proc.C)*upper_bound_factor, Nsteps)';
    results{pipe_loop_idx}.Dmodif.vals     = linspace(min(data_proc.D)*lower_bound_factor, max(data_proc.D)*upper_bound_factor, Nsteps)';
    results{pipe_loop_idx}.Sigmamodif.vals = linspace(min(data_proc.sigma)*lower_bound_factor,max(data_proc.sigma)*upper_bound_factor, Nsteps)';
   

    for IDX = 1:Nsteps
        % Run variation of A on pipe pipelist(pipe_loop_idx) % ================================
        [PRTgrv,PRTf] = run_simulation_func( ...
                            results{pipe_loop_idx}.Amodif.vals(IDX) ,...
                            data_proc.B(sample_select) ,...
                            data_proc.C(sample_select) ,...
                            data_proc.D(sample_select) ,...
                            data_proc.sigma(sample_select),...
                            Tend, ValveRampInit, ValveRampEnd, tvec);
        
            MX_PRTgrv(IDX,1)  = PRTgrv;
            MX_PRTf(IDX,1)    = PRTf;
        
        % Run variation of B        % ================================% ================================
        [PRTgrv,PRTf] = run_simulation_func(...
                        data_proc.Amax(sample_select),...
                        results{pipe_loop_idx}.Bmodif.vals(IDX),...
                        data_proc.C(sample_select),...
                        data_proc.D(sample_select),...
                        data_proc.sigma(sample_select) , ...
                        Tend, ValveRampInit, ValveRampEnd, tvec);

            MX_PRTgrv(IDX,2)  = PRTgrv;
            MX_PRTf(IDX,2)    = PRTf;

        
        % Run variation of C % ================================% ================================
        [PRTgrv,PRTf] = run_simulation_func(...
                data_proc.Amax(sample_select),...
                data_proc.B(sample_select),...
                results{pipe_loop_idx}.Cmodif.vals(IDX),...
                data_proc.D(sample_select), ...
                data_proc.sigma(sample_select) , ...
                Tend, ValveRampInit, ValveRampEnd, tvec);

            MX_PRTgrv(IDX,3)  = PRTgrv;
            MX_PRTf(IDX,3)    = PRTf;
            

        % Run variation of D % ================================% ================================
        [PRTgrv,PRTf] = run_simulation_func(...
                data_proc.Amax(sample_select),...
                data_proc.B(sample_select),...
                data_proc.C(sample_select), ...        
                results{pipe_loop_idx}.Dmodif.vals(IDX),...                
                data_proc.sigma(sample_select) , ...
                Tend, ValveRampInit, ValveRampEnd, tvec);

            MX_PRTgrv(IDX,4)  = PRTgrv;
            MX_PRTf(IDX,4)    = PRTf;
            

        % Run variation of Sigma % ================================% ================================
        [PRTgrv,PRTf] = run_simulation_func(...
                    data_proc.Amax(sample_select),...
                    data_proc.B(sample_select),...
                    data_proc.C(sample_select), ...
                    data_proc.D(sample_select),...                
                    results{pipe_loop_idx}.Sigmamodif.vals(IDX) ,...
                    Tend, ValveRampInit, ValveRampEnd, tvec);

            MX_PRTgrv(IDX,5)  = PRTgrv;
            MX_PRTf(IDX,5)    = PRTf;
            

    end 
    
    % Save ABCDsigma variational results for pipe pipelist(pipe_loop_idx)
    results{pipe_loop_idx}.results.MX_PRTgrv  = MX_PRTgrv;
    results{pipe_loop_idx}.results.MX_PRTf    = MX_PRTf;
    

end 
fprintf('Simulation took %1.3f\n',toc(tinit));


    % ===================================
    %           Plot results
    % ===================================

% ===============================================================

lgd = cell(length(pipelist), 1);
for idx = 1 : length(pipelist)
    lgd{idx} = sprintf('Sample %d', pipelist(idx));
end

% PARAMETER =====  A ========
figure(1); clf; 
%
axhA(1) = subplot(1,3,1); hold on;
for idx=1:length(pipelist)
    plot( 1e3*results{idx}.Amodif.vals,  1e3*results{idx}.results.MX_PRTgrv(:,1)  );
end
maxprtgrv      = [];
A_at_maxprtgrv = [];
for idx = 1:length(pipelist)-1
    [maxprtgrv(idx), idxtime] = max(results{idx}.results.MX_PRTgrv(:,1));
    A_at_maxprtgrv(idx)       = results{idx}.Amodif.vals(idxtime);
end
fitA = polyfit(A_at_maxprtgrv, maxprtgrv, 1);
param_range = [ results{1}.Amodif.vals(1) , results{1}.Amodif.vals(end) ];
plot(1e3*param_range, 1e3*polyval(fitA,param_range), '--k');
xlabel('A [ms]'); ylabel('[ms]'); title('PRTgrv'); grid on; box on;

for idx =1:length(pipelist) % Overlay the sample mean values
    plot(1e3*data_proc.Amax(pipelist(idx)), ...
        interp1( 1e3*results{idx}.Amodif.vals,  1e3*results{idx}.results.MX_PRTgrv(:,1), 1e3*data_proc.Amax(pipelist(idx) )),...
        'ok');
end


legend(lgd);
% -----------------------------
        axhA(2) = subplot(1,3,2); hold on;
        for idx=1:length(pipelist)
            plot( 1e3*results{idx}.Amodif.vals,  1e3*results{idx}.results.MX_PRTf(:,1)  );
        end
        maxprtf      = [];
        A_at_maxprtf = [];
        for idx = 1:length(pipelist)-1
            [maxprtf(idx), idxtime] = max(results{idx}.results.MX_PRTf(:,1));
            A_at_maxprtf(idx)       = results{idx}.Amodif.vals(idxtime);
        end
        fitA = polyfit(A_at_maxprtf, maxprtf, 1);        
        plot(1e3*param_range, 1e3*polyval(fitA,param_range), '--k');
        xlabel('A [ms]'); ylabel('[ms]'); title('PRTf');grid on; box on;
        for idx=1:length(pipelist)
            plot(1e3*data_proc.Amax(pipelist(idx)),...
                interp1(1e3*results{idx}.Amodif.vals,  1e3*results{idx}.results.MX_PRTf(:,1) ,1e3*data_proc.Amax(pipelist(idx))), ...
                'ok');
        end
        legend(lgd);
% -----------------------------
axhA(3) = subplot(1,3,3); hold on;
for idx=1:length(pipelist)
    try
    plot( 1e3*results{idx}.Amodif.vals,  results{idx}.results.MX_Pf_Pgrv(:,1)  );
    end
end
ylim([0 1]);
xlabel('A [ms]'); title('Pf/Pgrv Targ'); grid on; box on;

legend(lgd);
linkaxes(axhA,'x');
drawnow();

% PARAMETER =====  B ========
figure(2); clf; 
%
axhB(1) = subplot(1,3,1); hold on;
for idx=1:length(pipelist)
    plot( 1e3*results{idx}.Bmodif.vals,  1e3*results{idx}.results.MX_PRTgrv(:,2)  );
end
xlabel('B [ms]'); ylabel('[ms]'); title('PRTgrv'); 
grid on; box on;
for idx = 1 : length(pipelist)
    plot(1e3*data_proc.B(pipelist(idx)),...
        interp1(1e3*results{idx}.Bmodif.vals,  1e3*results{idx}.results.MX_PRTgrv(:,2) ,1e3*data_proc.B(pipelist(idx))),...
        'ok');
end
legend(lgd);
ylim(YLIMS);
% -----------------------------
        axhB(2) = subplot(1,3,2); hold on;
        for idx=1:length(pipelist)
            plot( 1e3*results{idx}.Bmodif.vals,  1e3*results{idx}.results.MX_PRTf(:,2)  );
        end
        ylim(YLIMS);
        xlabel('B [ms]'); ylabel('[ms]'); title('PRTf'); grid on; box on;
        for idx = 1 : length(pipelist)
            plot(1e3*data_proc.B(pipelist(idx)),...
                interp1(1e3*results{idx}.Bmodif.vals,  1e3*results{idx}.results.MX_PRTf(:,2),1e3*data_proc.B(pipelist(idx))),...
                'ok');
        end
        legend(lgd);
% -----------------------------
axhB(3) = subplot(1,3,3); hold on;
for idx=1:length(pipelist)
    try
    plot( 1e3*results{idx}.Bmodif.vals,  results{idx}.results.MX_Pf_Pgrv(:,2)  );
    end
end
ylim([0 1]);
xlabel('B [ms]'); title('Pf/Pgrv Targ'); grid on; box on; 
legend(lgd);
linkaxes(axhB,'x');
drawnow();



% PARAMETER =====  C ========
figure(3); clf; 
%
axhC(1) = subplot(1,3,1); hold on;
for idx=1:length(pipelist)
    plot( 1e3*results{idx}.Cmodif.vals,  1e3*results{idx}.results.MX_PRTgrv(:,3)  );
end
xlabel('C [ms]'); ylabel('[ms]'); title('PRTgrv'); grid on; box on; ylim(YLIMS);
for idx =1:length(pipelist)
    plot(1e3*data_proc.C(pipelist(idx)),...
        interp1(1e3*results{idx}.Cmodif.vals,  1e3*results{idx}.results.MX_PRTgrv(:,3),1e3*data_proc.C(pipelist(idx))),...
        'ok');
end
legend(lgd);
% -----------------------------
        axhC(2) = subplot(1,3,2); hold on;
        for idx=1:length(pipelist)
            plot( 1e3*results{idx}.Cmodif.vals,  1e3*results{idx}.results.MX_PRTf(:,3)  );
        end
        xlabel('C [ms]'); ylabel('[ms]'); title('PRTf'); grid on; box on;
        for idx=1:length(pipelist)
            plot(1e3*data_proc.C(pipelist(idx)),...
                interp1( 1e3*results{idx}.Cmodif.vals,  1e3*results{idx}.results.MX_PRTf(:,3),1e3*data_proc.C(pipelist(idx))),...
                'ok');
        end
        ylim(YLIMS);        
        legend(lgd);
% -----------------------------
axhC(3) = subplot(1,3,3); hold on;
for idx=1:length(pipelist)
    % plot( 1e3*1./(results{idx}.Cmodif.vals),  results{idx}.results.MX_Pf_Pgrv(:,3) )  );
    % plot( 1e3*1./(results{idx}.Cmodif.vals).^2,  1./(results{idx}.results.MX_Pf_Pgrv(:,3) )-1  );
    try
    plot( log10( results{idx}.Cmodif.vals ),  log10(  1./(results{idx}.results.MX_Pf_Pgrv(:,3)) -1  )  );    
    end
end
X = linspace(-3.13,-2.56);
Y = 2*X - mean(2*log10(data_proc.D))-0.1 ;
plot(X,Y,'--k');

% % % % plot([],+2*[],'-k'); % BUT with SmALL OFFSET in order to not hide the other curves 
xlabel('C [s] (log_{10})'); title('(Pgrv^o/Pf^o -1)'); grid on; box on;
ylabel('log_{10}');

lgd{length(pipelist)+1} = 'Slope +2x';
legend(lgd);
linkaxes(axhC([1,2]),'x');
drawnow();

% PARAMETER =====  D ========
figure(4); clf; 
%
axhD(1) = subplot(1,3,1); hold on;
for idx=1:length(pipelist)
    plot( 1e3*results{idx}.Dmodif.vals,  1e3*results{idx}.results.MX_PRTgrv(:,4)  );
end
xlabel('D [ms]'); ylabel('[ms]'); title('PRTgrv'); grid on; box on; ylim(YLIMS);
for idx=1:length(pipelist)
    plot(1e3*data_proc.D(pipelist(idx)),...
        interp1(1e3*results{idx}.Dmodif.vals,  1e3*results{idx}.results.MX_PRTgrv(:,4),1e3*data_proc.D(pipelist(idx))),...
        'ok');
end
legend(lgd);
% -----------------------------
        axhD(2) = subplot(1,3,2); hold on;
        for idx=1:length(pipelist)
            plot( 1e3*results{idx}.Dmodif.vals,  1e3*results{idx}.results.MX_PRTf(:,4)  );
        end        
        xlabel('D [ms]'); ylabel('[ms]'); title('PRTf'); grid on; box on;ylim(YLIMS);
        for idx=1:length(pipelist)
            plot(1e3*data_proc.D(pipelist(idx)),...
                interp1(1e3*results{idx}.Dmodif.vals,  1e3*results{idx}.results.MX_PRTf(:,4),1e3*data_proc.D(pipelist(idx))),...
                'ok');
        end
        legend(lgd);
% -----------------------------
axhD(3) = subplot(1,3,3); hold on;
for idx=1:length(pipelist)
    % plot( 1e3*1./(results{idx}.Dmodif.vals),  results{idx}.results.MX_Pf_Pgrv(:,4)  );
    % plot( 1e3*1./(results{idx}.Cmodif.vals).^2,  1./(results{idx}.results.MX_Pf_Pgrv(:,4) )-1  );
    try
    plot( log10( results{idx}.Dmodif.vals),  log10(  1./(results{idx}.results.MX_Pf_Pgrv(:,4)) -1 )  );
    end
end


X = linspace(-3.07,-2.7);
Y = -2*X + mean(2*log10(data_proc.C))-0.15 ;
plot(X,Y,'--k');
lgd{length(pipelist)+1} = 'Slope -2x';


xlabel('D [s] (log_{10})'); 
title('(Pgrv^o/Pf^o -1)'); grid on; box on;
ylabel('(log_{10})');
legend(lgd);
linkaxes(axhD([1,2]),'x');
drawnow();


% PARAMETER =====  Sigma ========
figure(5); clf; 
%
axhsig(1) = subplot(1,3,1); hold on;
for idx=1:length(pipelist)
    plot( results{idx}.Sigmamodif.vals,  1e3*results{idx}.results.MX_PRTgrv(:,5)  );
end
xlabel('\Sigma'); ylabel('[ms]'); title('PRTgrv'); grid on; box on;ylim(YLIMS);
legend(lgd);
% -----------------------------
        axhsig(2)=subplot(1,3,2); hold on;
        for idx=1:length(pipelist)
            plot( results{idx}.Sigmamodif.vals,  1e3*results{idx}.results.MX_PRTf(:,5)  );
        end        
        xlabel('\Sigma'); ylabel('[ms]'); title('PRTf'); grid on; box on; ylim(YLIMS);
        legend(lgd);
% -----------------------------
axhsig(3)=subplot(1,3,3); hold on;
for idx=1:length(pipelist)
    try
    plot( results{idx}.Sigmamodif.vals,  results{idx}.results.MX_Pf_Pgrv(:,5)  );
    end
end
ylim([0 1]);
xlabel('\Sigma'); title('Pf/Pgrv Targ'); grid on; box on; ylim([0 1]);
legend(lgd);
linkaxes(axhsig, 'x');
