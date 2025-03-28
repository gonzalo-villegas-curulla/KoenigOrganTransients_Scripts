clc, clear;
load data_proc.mat


% %%%%%%%%%%%%%%%%   CUSTOM USER PARAMETERS  %%%%%%%%%%%%%%%%%%%

pipelist = [8]; % \in[1,22] (DO NOT INCLUDE MORE THAN 3 PIPES at a time)

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

LOWER = 1e-3;
UPPER = 20e-3;
step  = 0.3e-3;

rang  = (LOWER : step : UPPER)';

% Parameter value allocation and simulation run:

tinit = tic;
results = cell(length(pipelist),1);
for pipe_loop_idx = 1:length(pipelist) % ===================================================================================================================

    sample_select = pipelist(pipe_loop_idx);
    
    results{pipe_loop_idx}.Amodif.vals     = rang;
    results{pipe_loop_idx}.Bmodif.vals     = rang;   

    for IDX = 1:length(results{pipe_loop_idx}.Amodif.vals)
        for JDX = 1 : length(results{pipe_loop_idx}.Bmodif.vals)
            [PRTgrv,PRTf,pf_over_pgrv_targ, flag_error] = run_simulation( ...
                                1/results{pipe_loop_idx}.Amodif.vals(IDX) ,...
                                1/results{pipe_loop_idx}.Bmodif.vals(JDX) ,...
                                1/data_proc.C(sample_select) ,...
                                1/data_proc.D(sample_select) ,...
                                data_proc.sigma(sample_select),...
                                Tend, ValveRampInit, ValveRampEnd, tvec);
            if flag_error
                MX_PRTgrv(IDX,JDX)  = nan;
                MX_PRTf(IDX,JDX)    = nan;
                MX_Pf_Pgrv(IDX,JDX) = nan;
            else
                MX_PRTgrv(IDX,JDX)  = PRTgrv;
                MX_PRTf(IDX,JDX)    = PRTf;
                MX_Pf_Pgrv(IDX,JDX) = pf_over_pgrv_targ;
            end
        end
        %
    end 
    % Save 
    results{pipe_loop_idx}.results.MX_PRTgrv  = MX_PRTgrv;
    results{pipe_loop_idx}.results.MX_PRTf    = MX_PRTf;
    results{pipe_loop_idx}.results.MX_Pf_Pgrv = MX_Pf_Pgrv;

end 
fprintf('Simulation took %1.3f\n',toc(tinit));


    % ===================================
    %           Plot results
    % ===================================

for idx = 1 : length(results{1}.results.MX_PRTgrv)
    MXtmp = zeros(size(results{1}.results.MX_PRTgrv));
    MXtmp(idx,:) = results{1}.results.MX_PRTgrv(idx,:);
    MXtmp(:,idx) = results{1}.results.MX_PRTgrv(:,idx);
    [aaa,bbb] = max(MXtmp,[],'all','omitnan');
    wheremaxgrv(idx) = bbb(end);
end
for idx = 1 : length(results{1}.results.MX_PRTf)
    MXtmp = zeros(size(results{1}.results.MX_PRTf));
    MXtmp(idx,:) = results{1}.results.MX_PRTf(idx,:);
    MXtmp(:,idx) = results{1}.results.MX_PRTf(:,idx);
    [aaa,bbb] = max(MXtmp,[],'all','omitnan');
    wheremaxf(idx) = bbb(end);
end
[AA,BB] = meshgrid(1e3*results{1}.Amodif.vals, 1e3*results{1}.Bmodif.vals);

% SLOPE PRTgrv max / A=B
mask_unique = 1:length(unique(BB(wheremaxgrv)));
slopemaxgrv = BB(wheremaxgrv(mask_unique))./AA(wheremaxgrv(mask_unique));

% SLOPE PRTf max / A=B
mask_unique = 1:length(unique(BB(wheremaxf)));
slopemaxf = BB(wheremaxf(mask_unique))./AA(wheremaxf(mask_unique));


% ===============================================================

lgd = cell(length(pipelist), 1);
for idx = 1 : length(pipelist)
    lgd{idx} = sprintf('Sample %d', pipelist(idx));
end


%% =====  PRT ========

RESOL = 10;

figure(1); clf; 
contourf(AA,BB, 1e3*results{1}.results.MX_PRTgrv, RESOL);
xlabel('A [ms]'); ylabel('B [ms]'); 
title(sprintf('PRTgrv [ms]'));%. Slope: %1.3f',median(slopemaxgrv)));
colorbar(); axis equal;
% hold on;
% plot(1e3*[LOWER,UPPER],1e3*[LOWER,UPPER], 'r', 'linewidth', 2);
% plot(AA(wheremaxgrv), BB(wheremaxgrv), 'om');


%% GVC KEEP [2025/03/17]
figure(2);clf;
plot( min(AA(:), BB(:)), ...
    1e3*results{1}.results.MX_PRTgrv(:) ,...
    '.');
hold on;
plot( min(AA(:), BB(:)) + 0.07, ...
    1e3*results{1}.results.MX_PRTf(:) ,...
    '.');
% axis equal;
grid on;
plot([0,20],[0,20],'-k', 'linewidth', 1.5);
xlabel('min(A,B) [ms]');
ylabel('PRT [ms]');


legend('PRT_{grv}','PRT_f');

%% [OK]

figure(3); clf;
plot(...
    BB(:),...
    results{1}.results.MX_PRTf(:)./results{1}.results.MX_PRTgrv(:),...
    '.');
grid on;
ylim([0 2]);
xlabel('B [ms]');
ylabel('PRT_f/PRT_{grv}');
%%

    figure();
    
    
    subplot(122);
    imagesc(1e3*results{1}.Amodif.vals, 1e3*results{1}.Bmodif.vals, 1e3*results{1}.results.MX_PRTf); ax=gca; ax.YDir = 'normal';
    xlabel('A [ms]'); ylabel('B [ms]'); 
    title(sprintf('PRTf [ms]. Slope: %1.3f',median(slopemaxf)));
    colorbar(); axis equal;
    hold on;
    plot(1e3*[LOWER,UPPER],1e3*[LOWER,UPPER], 'r', 'linewidth', 2);
    plot(AA(wheremaxf), BB(wheremaxf), 'om');
    ylim(1e3*[LOWER,UPPER]);
    figure(2); clf; 
    subplot(121);
    imagesc(1e3*results{1}.Amodif.vals, 1e3*results{1}.Bmodif.vals, 1e3*results{1}.results.MX_PRTgrv ); ax=gca; ax.YDir = 'normal';
    xlabel('A [ms]'); ylabel('B [ms]'); 
    title(sprintf('PRTgrv [ms]. Slope: %1.3f',median(slopemaxgrv)));
    colorbar(); axis equal;
    hold on;
    plot(1e3*[LOWER,UPPER],1e3*[LOWER,UPPER], 'r', 'linewidth', 2);
    plot(AA(wheremaxgrv), BB(wheremaxgrv), 'om');
    ylim(1e3*[LOWER,UPPER]);
    %
    subplot(122);
%%
figure();
    plot( min( results{1}.Amodif.vals(:), results{1}.Bmodif.vals(:) ) ,...
        results{1}.results.MX_PRTf(:), '.');
    hold on;
    % plot( min(1e3*results{1}.Amodif.vals,1e3*results{1}.Bmodif.vals),...
    %     results{1}.results.MX_PRTgrv(:), '.');
    legend();
%%


    % xlabel('A [ms]'); ylabel('B [ms]'); 
    % title(sprintf('PRTf [ms]. Slope: %1.3f',median(slopemaxf)));
    colorbar(); axis equal;
    hold on;
    plot(1e3*[LOWER,UPPER],1e3*[LOWER,UPPER], 'r', 'linewidth', 2);
    plot(AA(wheremaxf), BB(wheremaxf), 'om');
    ylim(1e3*[LOWER,UPPER]);



% ================================== PRT vs min(A,B)
YLIMS = 1.05*[0, 1e3*max( max(results{1}.results.MX_PRTgrv,[],'all'),max(results{1}.results.MX_PRTf,[],'all') ) ];
LW = 1.5;



%% [NO]
figure();
plot( BB(:), ...
    exp(results{1}.results.MX_PRTf(:)./results{1}.results.MX_PRTgrv(:)) ,...
    '.');
%%



%%
figure(10);clf;
subplot(2,2,1);
plot( min(AA(:), BB(:)), ...
    1e3*results{1}.results.MX_PRTgrv(:) ,...
    'o');
grid on;
xlabel('min(A,B) [ms]');
ylabel('PRT [ms]');
axis equal;
hold on;
plot(1e3*[LOWER,UPPER],1e3*[LOWER,UPPER],'-r', 'linewidth',LW);
title('Groove');
ylim(YLIMS);
%
subplot(2,2,2);
plot( max(AA(:), BB(:)), ...
    1e3*results{1}.results.MX_PRTgrv(:) ,...
    'o');
grid on;
xlabel('max(A,B) [ms]');
ylabel('PRT [ms]');
axis equal;
hold on;
plot(1e3*[LOWER,UPPER],1e3*[LOWER,UPPER],'-r', 'linewidth',LW);
title('Groove');
ylim(YLIMS);
%
subplot(2,2,3);
plot( min(AA(:),BB(:)),...
    1e3*results{1}.results.MX_PRTf(:),...
    'o');
xlabel('min(A,B) [ms]');
axis equal;
hold on;
plot(1e3*[LOWER,UPPER],1e3*[LOWER,UPPER],'-r', 'linewidth',LW);
title('Foot'); grid on;
ylabel('PRT [ms]');
ylim(YLIMS);
%
subplot(2,2,4);
plot(max(AA(:),BB(:)), ...
    1e3*results{1}.results.MX_PRTf(:),...
    'o');
axis equal;
hold on;
plot(1e3*[LOWER,UPPER], 1e3*[LOWER,UPPER], '-r', 'linewidth',LW);
title('Foot'); grid on;
ylabel('PRT [ms]');
ylim(YLIMS);


% ========
% figure(11); clf;
% surf( 1e3*AA, 1e3*BB, 1e3*results{1}.results.MX_PRTgrv, 'linestyle',  'none');

%%


figure();

plot(results{1}.results.MX_PRTgrv(:), ...
    results{1}.results.MX_PRTf(:));
%%
figure();
plot( results{1}.results.MX_PRTgrv(:),  ...
    (results{1}.results.MX_PRTgrv(:)-results{1}.results.MX_PRTf(:)  ) , 'o');

%%
figure();

plot( results{1}.results.MX_PRTf(:),  ...
    max(results{1}.results.MX_PRTf(:)./results{1}.results.MX_PRTgrv(:) ),...
    'o');
grid on;
xlabel();
ylabel('PRTf/PRTgrv');



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
                figure(10);clf;
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
