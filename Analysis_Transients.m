try
    cd /run/media/gvc/ExtremeSSD/OrganPipe2023-2024/DataTransients/
end 
clc; clear; 
fs = 51.2e3;
dt = 1/fs;
    videoObj = VideoWriter('output_video.avi'); 
    videoObj.FrameRate = 12; 
    open(videoObj);

files=dir('A*.mat');
for BIGIDX = 5 :  length(files)  % ==== MAIN LOOP ====

clc; clearvars -except BIGIDX files videoObj
% close all;
tinit = tic();
fs    = 51.2e3;
dt    = 1/fs;

filename = files(BIGIDX).name;
thedata  = load(filename);
tvec     = thedata.tvec;

x = thedata.PR_data.pressureData;
for idx = 1 : 6
    x(idx,:) = x(idx,:)-mean(x(idx,1:100)); % Correct AUTO-ZERO detection of contidioners
end

% Retrieve key movement for reference indexing
[KeyDownIdx,KeyUpIdx,KeyMovingTime,DurNotesInS, VelPeakIdxPos,VelPeakIdxNeg] = DetectVelocityPeaks_func(x(6,:),tvec,filename);

if length(VelPeakIdxPos)~=length(VelPeakIdxNeg)
    VelPeakIdxPos = VelPeakIdxNeg + fix(fs*3);
end
NT = length(KeyDownIdx);


% INITS/ VECTOR ALLOCATIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pgroove_targ  = zeros(NT,1);
Pfoot_targ    = zeros(NT,1);
Ppipe_targ    = zeros(NT,1);
PpalletB_targ = zeros(NT,1);

t20groove = zeros(NT,1);
t20foot   = zeros(NT,1);
t20mouth  = zeros(NT,1);
Area1     = zeros(NT,1);
Area2     = zeros(NT,1);

afit      = zeros(NT,1);
bfit      = zeros(NT,1);
cfit      = zeros(NT,1);
dfit      = zeros(NT,1);
gofr2     = zeros(NT,1);

A2max_over_A1simult = zeros(NT,1);
A2max_over_A1target = zeros(NT,1);
A2max_over_A2target = zeros(NT,1);
a2max_vec = zeros(NT,1);
pf_at_a2max = zeros(NT,1);
pm_at_a2max = zeros(NT,1);
max_a2_over_a1 = zeros(NT,1);


% Retrieve transient location running (DetectVelocityPeaks.m) ============
lk_foot     = KeyDownIdx;
lk_foot_end = lk_foot' + fix(fs*DurNotesInS);

    if 0 % PLOT signals and onsets???
        plotquick; hold on; plot(lk_foot,zeros(length(lk_foot)),'ob');
        plot(lk_foot_end, zeros(length(lk_foot_end)),'*r');
    end


% A prediction of transient length to allocate
LENretro = fix(0.100*fs); 
LENpost  = fix(2.2*fs + 0*0.600*fs);  % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< PUT IT BACK TO 0.600f*s

AvTimeDur     = 0.500; % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segment and loop over all transients of current file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fighand = figure(20);clf; 
axhand = axes(fighand); cla;

for idx = 1 : length(lk_foot)
        
        % STORE TRANSIENTS =================================
        
        try
            midpoint_idx  = fix( 0.5*(lk_foot(idx) + lk_foot_end(idx)) );
            if ~isnan(midpoint_idx)
                midpoint_span = midpoint_idx:midpoint_idx + fix(fs*AvTimeDur);
            else
                midpoint_span = lk_foot(idx)+fix(fs*1.5):lk_foot(idx)+fix(fs*1.5) + fix(fs*AvTimeDur);
                flag = true;
            end
        catch
            midpoint_span = lk_foot(idx):lk_foot(idx)+fix(fs*1.5) + fix(fs*AvTimeDur);
        end

        
        % Both masks start at the same point, otherwise, idxmess for t20%
        mask      = lk_foot(idx)-LENretro : lk_foot(idx)+LENpost;
        if ~flag
            mask_pipe = lk_foot(idx)-LENretro : lk_foot_end(idx)+LENpost;
        else 
            mask_pipe = lk_foot(idx)-LENretro : lk_foot(idx)+fix(fs*3)+LENpost;
        end
        
        
        carefulmask = lk_foot(idx) : lk_foot(idx)+fix(3*fs);
        ll = fix(length(carefulmask)*0.5): fix(length(carefulmask)*0.75);
        carefulmask = carefulmask( ll );
        
        
        % PALLET BOX  /!\  ~~~~~~~~~~~~~~~~~
        PpalletB_targ(idx) = mean(x(2,carefulmask));
        
        
        % GROOVE /!\  ~~~~~~~~~~~~~~~~~
        
        tmp_gr = x(3, mask );
        OFFSET = mean(tmp_gr(1:100));
        tmp_gr = tmp_gr - OFFSET;
        groove_trans{idx} = tmp_gr;
        
      
        Pgroove_targ(idx) = mean(x(3, carefulmask ))-OFFSET;        
        t20idx            = find(tmp_gr/Pgroove_targ(idx) >0.2, 1, 'first');
        t20idxgr = t20idx;
        t80idx            = find(tmp_gr/Pgroove_targ(idx) >0.8, 1, 'first');        
        % t20groove(idx) = dt*(t20idx-LENretro -lk_foot(idx)); % in [s]
        t20groove(idx)    = dt*(t20idx-LENretro ); % in [s]
        PRTgroove(idx)    = dt*(t80idx-t20idx);
        
        if 0
            figure(12);clf;plot(tmp_gr); hold on;
            plot([1,length(tmp_gr)],[1,1]*Pgroove_targ(idx)*0.8,'--k');
            plot([1,length(tmp_gr)],[1,1]*Pgroove_targ(idx)*0.2,'--k');
            pause();
        end                        
        

        % FOOT  /!\  ~~~~~~~~~~~~~~~~~
        tmp_ft          = x(4, mask );
        OFFSET          = mean(tmp_ft(1:100));
        tmp_ft          = tmp_ft - OFFSET;
        foot_trans{idx} = tmp_ft;
        
        Pfoot_targ(idx) = mean(x(4, carefulmask  ))- OFFSET;
        t20idx = find(tmp_ft/Pfoot_targ(idx) >0.2, 1, 'first');
        t80idx = find(tmp_ft/Pfoot_targ(idx) >0.8, 1, 'first');
        t20foot(idx) = dt*(t20idx-LENretro  ); % in seconds
        PRTfoot(idx) = dt*(t80idx  - t20idx );
        
                % RADIATED SOUND/PIPE/MOUTH  /!\  ~~~~~~~~~~~~~~~~~
        
        tmp_wholepipe    = x(5, mask_pipe);
        tmp_wholepipe    = tmp_wholepipe - mean(tmp_wholepipe);
        pipedatavec{idx} = tmp_wholepipe ;
            
        tmp_rad        = x(5, mask);
        tmp_rad        = tmp_rad - mean(tmp_rad(1:100));
        pipetrans{idx} = tmp_rad;
        
        
        % F1 and T1
        SSmask = lk_foot(idx) : lk_foot(idx)+fix(3*fs); % until 3 seconds after key-down
        ll     = fix(length(SSmask)*0.5): fix(length(SSmask)*0.75);
        SSmask = SSmask( ll );
        midpipedata = x(5, SSmask )'; midpipedata = midpipedata(:);

        samplepipe = filename(2:3);

        switch samplepipe
            case '04'
                % fprintf("f1 YIN 155Hz, ")
                KnownFreq = 155.56; 
                P.sr      = fs;
                P.minf0 = 0.9*KnownFreq;
                P.maxf0 = 1.1*KnownFreq;
                Ryin    = yin( midpipedata  ,P);
                f1estim = mean(Ryin.f0_Hz,'omitnan');
                T1estim = mean(1./f1estim,'omitnan');
            case '05'
                % fprintf("f1 YIN 164Hz, ")
                KnownFreq = 164.81; 
                P.sr      = fs;
                P.minf0 = 0.9*KnownFreq;
                P.maxf0 = 1.1*KnownFreq;
                Ryin    = yin( midpipedata  ,P);
                f1estim = mean(Ryin.f0_Hz,'omitnan');
                T1estim = mean(1./f1estim,'omitnan');

            case '17'
                % fprintf("f1 YIN 329Hz, ")
                KnownFreq = 329.0;
                P.sr      = fs;
                P.minf0 = 0.9*KnownFreq;
                P.maxf0 = 1.1*KnownFreq;
                Ryin    = yin( midpipedata  ,P);
                f1estim = mean(Ryin.f0_Hz,'omitnan');
                T1estim = mean(1./f1estim,'omitnan');

            case '24'
                % fprintf("f1 YIN 493Hz, ")
                KnownFreq = 493.0;
                P.sr      = fs;
                P.minf0 = 0.9*KnownFreq;
                P.maxf0 = 1.1*KnownFreq;
                Ryin    = yin( midpipedata  ,P);
                f1estim = mean(Ryin.f0_Hz,'omitnan');
                T1estim = mean(1./f1estim,'omitnan');
                
            otherwise
                % fprintf("f1 xcorr, ");
                [RR,lags] = xcorr(midpipedata);
                [~,locs]  = findpeaks(RR, 'MinPeakProminence',(max(RR)-std(RR)));
                T1estim   = mean(diff(locs))*dt;
                f1estim   = 1/T1estim;
        end

        f1(idx) = f1estim;
        TN = fix(T1estim*fs);

        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FOURIER DECOMPOSITION of Prad
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        try
            fundam = mean(f1,'omitnan');
        catch
            error('Could not establish Fundam [Hz] previously. What should I use?');
        end
        
        pipecurr     = pipedatavec{idx};
        [comp0,COMP] = DecomposeFourier_func(pipecurr, fs, fundam);
        COMP         = COMP'; % Fourier first pressure components

        run myenvelopedetection.m
        

        Ppipe_targ(idx) = TOTAL_target_amp; % Mouth radiated target pressure
        t20idxp         = find( enveltot/Ppipe_targ(idx)<0.2,1,'last');
        t80idxp         = find( enveltot/Ppipe_targ(idx)>0.8,1,'first');

        % t20mouth(idx)   = dt*(t20idxp-LENretro -lk_foot(idx));
        t20mouth(idx) = dt*(t20idxp-LENretro );
        PRTpipe(idx)  = dt*(t80idxp-t20idxp);  

        
        % Normalize by their steady-state values
        Fund     = envel_first/mean(envel_first(mask3));
        Second   = envel_second/mean(envel_second(mask3));
        Third    = envel_third/mean(envel_third(mask3));
        ll = t20idxp:  t80idxp; 
        try % Area integration of F2/F1@trans and F3/F1@transient
            Area1(idx) = trapz(dt, Second(ll)-Fund(ll)  ); % equally spaced in time (simplify to dt)
            Area2(idx) = trapz(dt, Third(ll)-Fund(ll)  );
        end


        % We try to find point-wise descriptors to compare a2 with a1:
        [a2max_val,a2max_idx]   = max(envel_second(1:length(tmp_ft)));
        a2targ = mean(envel_second(end-fix(50*fs/f1estim): end));
        A2max_over_A1simult(idx) = envel_second(a2max_idx)/envel_first(a2max_idx);
        A2max_over_A1target(idx) = a2max_val/mean(envel_first(end-fix(fs*0.500):end),'omitnan');
        A2max_over_A2target(idx) = a2max_val/target_second; 

        % For transient characterisation after fig6, CFA-Ernoult2016:
        pf_at_a2max(idx) = tmp_ft(a2max_idx);
        pm_at_a2max(idx) = pipecurr(a2max_idx);
        a2max_vec(idx) = a2max_val;
    
        % ======================================================
        % Compute and plot the ratio between a2 and a1 to find the max
        % between t20_foot and t80_foot, after having smoothed them with
        % butterworth() filtfilt() with a length of N periods (Def. 5).

        [bbb,aaa]=butter(4, (f1estim/5)/(fs/2));
        
        dat1 = filtfilt(bbb,aaa,envel_first);
        dat2 = filtfilt(bbb,aaa,envel_second);
        ratio_dats = dat2./dat1;

        max_a2_over_a1(idx) = max(ratio_dats(  t20idxp :  t80idxp+fix(50*PRTfoot(idx)*fs)   ));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %        VISUALIZATIONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if 0
            %%% PHASE-SPACE TRANSIENT Pgroove-Pfoot
            maskvid = 1 +0*t20idx:t80idx+fix(1*fs*PRTfoot(idx)) ;
                midmask = fix(0.5*length(tmp_ft)):fix(0.5*length(tmp_ft)) + fix(20*fs/f1estim);
            % plot(axhand, tmp_gr(maskvid)/Pgroove_targ(idx) ,   tmp_ft(maskvid)/Pfoot_targ(idx), 'linewidth',2 );                 
            plot(axhand, tmp_gr(midmask)/Pgroove_targ(idx) ,   tmp_ft(midmask)/Pfoot_targ(idx), 'linewidth',0.5 ); 
            plot(axhand, (tmp_ft(midmask)-Pfoot_targ(idx))/Ppipe_targ(idx),  tmp_rad(midmask)/Ppipe_targ(idx), 'linewidth',0.5);

            grid on; box on; hold(axhand, 'on');
                plot(axhand, tmp_gr(t20idxgr)/Pgroove_targ(idx)*[1,1],[0,0.5],'--r');
                plot(axhand, [0,1.0],tmp_ft(t20idx)/Pfoot_targ(idx)*[1,1],'--r')
                plot([1,1],[0,1],'--k'); plot([0,1],[1,1],'--k'); plot([0,1],[0,1],'--k');
            hold(axhand, 'off');

            xlabel('Pgroove'); ylabel('Pfoot'); legend('Pf/Pgr','t20','location','NorthWest');
            axis equal;xlim([0.7 1.3]); ylim([0.7 1.3]);
            title(sprintf([files(BIGIDX).name, ', trans num: ', num2str(idx)]), 'interpreter','none');        

            drawnow();
            frame = getframe(gcf);
            writeVideo(videoObj, frame);
        end

        if 1
            N1period = fix(fs/f1estim);
            HA       = fix(0.2*N1period);
            sliding_mask = [ 1 : 10*N1period] + 0*HA;
            while sliding_mask(end) < length(tmp_rad)
                xwin = (tmp_ft(sliding_mask)-Pfoot_targ(idx))/Ppipe_targ(idx);
                ywin = tmp_rad(sliding_mask)/Ppipe_targ(idx);
                if 1
                    plot(axhand, xwin, ywin, 'linewidth',0.5); hold on;
                    xlim(3.0*[-1,1]); ylim(3.0*[-1,1]);
                    plot([0,0],[-1,1],'--r');plot([-1,1],[0,0],'--r');plot([-1,1],[-1,1],'--r');plot([-1,1],[1,-1],'--r'); 
                else
                    plot3(axhand, xwin,[1:length(sliding_mask)]*dt ,ywin,'linewidth',1); hold on; 
                    % axis equal;
                    xlim(3.0*[-1,1]); zlim(3.0*[-1,1]);
                    view(141.42,21.9);
                end

                box on; grid on; 
                xlabel('(Pfoot-Pfoottarget) / Ppipetarget'); ylabel('Ppipe / Ppipetarget');
                title(sprintf([files(BIGIDX).name, ', trans num: ', num2str(idx)]), 'interpreter','none');        
                hold off;
                sliding_mask = sliding_mask + HA;
                
                drawnow();
                % pause();
            end
            
        end

        % ===================================================
        if 0
            figure(1);clf;hold on;

            plot([envel_first;envel_second]');grid on;
            % xlim([0.5,1.5]*1e4);
            xlim([5e3, t80idxp+fix(50*PRTfoot(idx)*fs)]);

            plot(dat1,'--k');
            plot(dat2,'--k');
            
            plot(t20idxp*[1,1],[0 a2max_val],'--k');plot(t80idxp*[1,1],[0 a2max_val],'--k');
            
            yyaxis right;plot(dat2./dat1,'g');ylim([0,1]*20); box on;
            title(sprintf([files(BIGIDX).name, ', trans num: ', num2str(idx)]), 'interpreter','none');              
            drawnow();
            % pause(1);
            frame = getframe(gcf);
            writeVideo(videoObj, frame);
        end
        % ===================================================

        if 0 % <<Plot envelopes, t20, t80, 1-2-3 harmonics, and total rad pressure>>

            TLIMS = [0 0.500];
            figure(23); clf;
            subplot(2,2,1);
            plot(time,pipecurr(1:hL) );
            hold on;
            plot(time, enveltot, 'r'); ylabel('Pressure MouthRad total');
            plot(time, envel_sum, 'g');
            plot([0,time(end)],Ppipe_targ(idx)*[1,1],'--k');
            plot([0,time(end)],0.8*Ppipe_targ(idx)*[1,1],'--k');
            plot([0,time(end)],0.2*Ppipe_targ(idx)*[1,1],'--k');
            plot([1,1]*time(t20idxp),[0,1.1*max(pipecurr) ],'--k');
            plot([1,1]*time(t80idxp),[0,1.1*max(pipecurr) ],'--k');
            xlim(TLIMS);
            rawLim = 1.1*max(pipecurr);
            ylim([-5 rawLim]);
            title(sprintf([files(BIGIDX).name, ', trans num: ', num2str(idx)]), 'interpreter','none');        
            
            subplot(2,2,2);
            plot(time, COMP(1,1:hL)); hold on; plot(time, envel_first, 'r');ylabel('Fundam'); ylim([-5,rawLim]);xlim([0 0.500]);
            
            subplot(2,2,3);
            plot(time,COMP(2,1:hL));hold on; plot(time, envel_second, 'r');ylabel('2nd Harmonic');xlabel('Time [s]');ylim([-5,rawLim]);xlim([0 0.500]);
        
            subplot(2,2,4);
            plot(time, COMP(3,1:hL));
            hold on; plot(time, envel_third, 'r');ylabel('3rd Harmonic');ylim([-5,rawLim]);xlim([0 0.500]);
            drawnow;
            % pause(2);            
        end

                
end % Close loop over all transients of current file
fprintf("\n");
% Period (found by yin sometimes):
T1 = 1./f1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FOOT PRESSURE FIT beta and nu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PREPARE FIT OPTIONS ===========================

ft = fittype( 'a./(1+d*exp(-b*(x-c))).^(1/d)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' , 'TolFun', 1e-12);
opts.Display = 'Off';

%                   (a)             (b)   (c)    (d)
%                  Ptarg            beta  t_o    nu
opts.Lower      = [200              300   1e-3   1e-4 ];
opts.Upper      = [600              10e3  2      10    ]; 
opts.StartPoint = [mean(Pfoot_targ) 300   0.01   0.01];
opts.Robust     = 'Bisquare'; % LAR, Off, Bisquare

% /////////////////////////////////////////////////////////


% ALIGN AND GROUP-PLOT  ===========================
PRECUT = fix(0.020*fs); % (Def. 0.015 s) Shift Data away from time zero
fprintf("Starting beta nu fits...");
fighand = figure(20);clf; 
axhand = axes(fighand); cla;

NumPRTs = 1; %length of data after t80foot for the fit

if 0 % for jdx = 1 : length(foot_trans) % LOOP OVER ALL TRANSIENTS of current file <<<<<<<<<<<<<<<<<<<<

    Ptar = Pfoot_targ(jdx);
    
    %                           (b)   (c)    (d)
    %                           beta  t_o    nu
    opts.Lower      = [0.6*Ptar 300   1e-3   3e-3 ];
    opts.Upper      = [2.0*Ptar 10e3  2      10    ]; 
    opts.StartPoint = [Ptar     300   0.01   0.5];
    opts.Robust     = 'Bisquare'; % LAR, Off, Bisquare

    init_idx = find( foot_trans{jdx}/Pfoot_targ(jdx)< 0.2, 1, 'last') - fix(fs*5e-3);
    end_idx  = find( foot_trans{jdx}/Pfoot_targ(jdx)>0.8, 1, 'first') + fix(NumPRTs*PRTfoot(jdx)*fs);

    % Fit model to data.
    yData = foot_trans{jdx}(init_idx : end_idx);
    xData = [0:length(yData)-1]*dt;
    [FitRes{jdx}, gof{jdx}] = fit( xData', yData', ft, opts );


    % Plot fit with data.
    if 0
        
        plot(axhand, FitRes{jdx}, xData, yData);
        legend('Pressure vs. Time', 'Fitted model', 'Location', 'SouthEast', 'Interpreter', 'none' );
        % Label axes
        xlabel( 'Time [s]', 'Interpreter', 'none' );
        ylabel( 'Foot pressure [Pa]', 'Interpreter', 'none' );
        grid on;
        
        title(sprintf([files(BIGIDX).name, ', trans num: ', num2str(jdx)]), 'interpreter','none');        
        drawnow;
        frame = getframe(gcf);
        writeVideo(videoObj, frame);
        % pause(0.5);
    end

    
    if 0 
    fprintf('Pfoot_targ %1.1f, Beta %1.1f, Nu %1.5f\n',FitRes{jdx}.a,FitRes{jdx}.b,FitRes{jdx}.d);
    end

    % Store fit data
    afit(jdx)  = FitRes{jdx}.a; % notPtarg
    bfit(jdx)  = FitRes{jdx}.b; % beta
    cfit(jdx)  = FitRes{jdx}.c;
    dfit(jdx)  = FitRes{jdx}.d; % nu
    gofr2(jdx) = gof{jdx}.rsquare;

    
end % CLOSE LOOP over all transients of current file
fprintf(" (Done)\n");

Ptargfit = afit;
betafit  = bfit;
delayfit = cfit;
nufit    = dfit;


% Reduced Jet Velocity in S-S if Wm is known ==========================

RJV = sqrt(2*Pfoot_targ(:)/1.2)./(thedata.PR_params.Wm*f1(:));
Wm  = thedata.PR_params.Wm;

% %%%%%% SAVE processed data and results %%%%%%%%%%

datafilename = filename(1:end-4);

if 0
    save(['./processed/' datafilename '_PROCESSED.mat'],'f1','betafit','nufit','Ptargfit','delayfit',...
        'Area1','Area2','RJV','Wm','fs','PpalletB_targ','Pgroove_targ','t20groove',...
        'PRTgroove','Pfoot_targ','t20foot','PRTfoot','Ppipe_targ','t20mouth','PRTpipe',...
        'KeyMovingTime','gofr2','A2max_over_A1simult','A2max_over_A1target',...
        'A2max_over_A2target','pf_at_a2max','pm_at_a2max','a2max_vec','max_a2_over_a1');    

    
    fprintf('Ellapsed time with this file: %1.2f \n', toc(tinit));
end


end % BIGIDX of all files opened

close(videoObj);