function [PRTg, PRTf] = run_simulation_func(...
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
    
    % [t_ode,y] = ode45(@(t_ode,y) solverA(t_ode, y,...
    [t_ode,y] = ode45(@(t_ode,y) solverOrgan(t_ode, y,...
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
    PRTg   = t90grv-t10grv;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 ACHTUNG:                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if pgrv(end)<pf(end)
        PRTg = NaN;
    end


    %
    t10f = tvec(find(pf/pf(end)<0.1,1,'last'));
    t90f = tvec(find(pf/pf(end)>0.9,1,'first'));
    PRTf = t90f-t10f;

    % fprintf('%1.2f %1.2f; \n ',1e3*PRTg, 1e3*PRTf);

    if 0 % Plot all time integrations

        figure(10);clf;
        LW = 1.5;
        plot(t_ode*1e3, y(:,1),'linewidth',LW);
        hold on;
        plot(t_ode*1e3, y(:,2),'linewidth',LW);                              
        % fprintf(sprintf('Max pgrv: %1.3f. Max pf: %1.3f. Nan grv %d. Nan f %d.\n',...
            % max(pgrv), max(pf), isnan(pgrv(end)),isnan(pf(end)) ));        
        grid on;
        xlabel('Time [ms]');ylabel('pressure [n.u.]'); 
        title(sprintf('Solution in time. Ta: %1.2f (ms); SIG: %1.1f',PASS_Amax,PASS_sigma_full));
        xlim([50 150]); ylim([0 1]);
        drawnow();
        % writeVideo(vobj, getframe(gcf) );                
        pause();

    end
end
%
function dydt = solverA(t_ode, y,A,B,C,D,sigMa_full, ValveRampInit, ValveRampEnd)

    omeg = omega_func(t_ode, ValveRampInit, ValveRampEnd);    
    dydt = zeros(2,1);
    dydt(1) = omeg*real(sqrt(2*(1-y(1))))/A   -   real(sqrt(y(1)*(2-2*omeg.^2*sigMa_full.^2) -2*y(2) + 2*omeg.^2*sigMa_full^2 ))/B;
    dydt(2) =      real(sqrt(y(1)*(2-2*omeg.^2*sigMa_full.^2) -2*y(2)+2*omeg.^2*sigMa_full.^2 ))/C  -  real(sqrt(2*y(2)))/D;

end