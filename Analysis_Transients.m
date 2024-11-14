clc; clear; 
fs = 51.2e3;
dt = 1/fs;
    videoObj = VideoWriter('output_video.avi'); 
    videoObj.FrameRate = 12; 
    open(videoObj);

files=dir('A*.mat');
for BIGIDX = 1 : length(files)  % ==== MAIN LOOP ====

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
dfit      = zeros(NT,1);
gofr2     = zeros(NT,1);

A2max_over_A1simult = zeros(NT,1);
A2max_over_A1target = zeros(NT,1);


% Retrieve transient location running (DetectVelocityPeaks.m) ============
lk_foot     = KeyDownIdx;
lk_foot_end = lk_foot' + fix(fs*DurNotesInS);

    if 0 % PLOT signals and onsets???
        plotquick; hold on; plot(lk_foot,zeros(length(lk_foot)),'ob');
        plot(lk_foot_end, zeros(length(lk_foot_end)),'*r');
    end


% A prediction of transient length to allocate
LENretro = fix(0.100*fs); 
LENpost  = fix(0.600*fs); 

AvTimeDur     = 0.500; % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segment and loop over all transients of current file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        tmp_ft = x(4, mask );
        OFFSET = mean(tmp_ft(1:100));
        tmp_ft = tmp_ft - OFFSET;
        foot_trans{idx} = tmp_ft;
        
        Pfoot_targ(idx) = mean(x(4, carefulmask  ))- OFFSET;
        t20idx = find(tmp_ft/Pfoot_targ(idx) >0.2, 1, 'first');
        t80idx = find(tmp_ft/Pfoot_targ(idx) >0.8, 1, 'first');
%         t20foot(idx) = dt*(t20idx-LENretro -lk_foot(idx)  ); % in seconds
        t20foot(idx) = dt*(t20idx-LENretro  ); % in seconds
        PRTfoot(idx) = dt*(t80idx  - t20idx );
        
        if 0
            figure(12);clf;plot(tmp_ft);  hold on;
            plot([1,length(tmp_ft)],[1,1]*Pfoot_targ(idx)*0.8,'--k');
            plot([1,length(tmp_ft)],[1,1]*Pfoot_targ(idx)*0.2,'--k');
        end

      
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
        if 1 % <<<<<< % A04 (sample2) & A05 (samp.3) fail: they should be 155.56Hz and 164.81
            [RR,lags] = xcorr(midpipedata);
            [~,locs]  = findpeaks(RR, 'MinPeakProminence',(max(RR)-std(RR)));
            T1estim   = mean(diff(locs))*dt;
            f1estim   = 1/T1estim;
        else
            fprintf("Xcorr failed for f1estimation, trying with YIN...\n")
            KnownFreq = 164.81; % A04 (sample2) & A05 (samp.3) fail: they should be 155.56Hz and 164.81
            P.sr      = fs;
            P.minf0 = 0.9*KnownFreq;
            P.maxf0 = 1.1*KnownFreq;
            Ryin    = yin( midpipedata  ,P);
            f1estim = mean(Ryin.f0_Hz,'omitnan');
            T1estim = mean(1./f1estim,'omitnan');
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
            Area1(idx) = trapz([0:length(ll)-1]*dt, Second(ll)-Fund(ll)  );
            Area2(idx) = trapz([0:length(ll)-1]*dt, Third(ll)-Fund(ll)  );
        end
        [a2max_val,a2max_idx] = max(envel_second);
        A2max_over_A1simult(idx) = envel_second(a2max_idx)/envel_first(a2max_idx);
        A2max_over_A1target(idx) = a2max_val/mean(envel_first(end-fix(fs*0.500):end),'omitnan');

        if 1 % <<Plot envelopes, t20, t80, 1-2-3 harmonics, and total rad pressure>>

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

% Period (found by yin sometimes):
T1 = 1./f1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FOOT PRESSURE beta and nu fit
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

 % lkjdfsgoiuhvtpo4q3qut34qiubytvel
 % bternklrwlukvtewku348io34i7yn3qq5qbvihubtwitrw4ibltr4wtviw7478w45ibyw4a5liw45iw4854wy5i4wyi5w45hlibw4va5lnhw5wal8i5aw4v5

for jdx = 1 : length(foot_trans) % LOOP OVER ALL TRANSIENTS of current file

    Ptar = Pfoot_targ(jdx);
    
    ft = fittype( 'Ptar./(1+d*exp(-b*(x-c))).^(1/d)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' , 'TolFun', 1e-12);
    opts.Display = 'Off';
    %                  (b)   (c)    (d)
    %                  beta  t_o    nu
    % opts.Lower      = [Ptar 300   1e-3   1e-4 ];
    % opts.Upper      = [Ptar 10e3  2      10    ]; 
    % opts.StartPoint = [Ptar 300   0.01   0.01];
    % opts.Robust     = 'Bisquare'; % LAR, Off, Bisquare

    %                  (b)   (c)    (d)
    %                  beta  t_o    nu
    opts.Lower      = [Ptar 300   1e-3   1. ];
    opts.Upper      = [Ptar 10e3  2      1. ]; 
    opts.StartPoint = [Ptar 300   0.01   1.];
    opts.Robust     = 'Bisquare'; % LAR, Off, Bisquare



    
    shiftonset       = find(foot_trans{jdx}/Pfoot_targ(jdx)>0.1,1,'first') ;

    footdatavec{jdx} = circshift(foot_trans{jdx} , -shiftonset + PRECUT ); 
    timevec{jdx}     = ([0:length(footdatavec{jdx})-1] - PRECUT)*dt  ;           

    TRstartidx = 1;
    TRendidx   = length(footdatavec{jdx})-shiftonset+PRECUT;

    % Fit model to data.
    % try
        xData = timevec{jdx}(1:TRendidx)';
        yData = footdatavec{jdx}(1:TRendidx);
        [FitRes{jdx}, gof{jdx}] = fit( xData, yData', ft, opts );
    % end

    % Plot fit with data.
    if 1
        figure(20);clf; 
        h = plot( FitRes{jdx}, xData, yData );
        legend( h, 'Pressure vs. Time', 'Fitted model', 'Location', 'NorthEast', 'Interpreter', 'none' );
        % Label axes
        xlabel( 'Time [s]', 'Interpreter', 'none' );
        ylabel( 'Foot pressure [Pa]', 'Interpreter', 'none' );
        grid on;
        
        title(sprintf([files(BIGIDX).name, ', trans num: ', num2str(jdx)]), 'interpreter','none');        
        % ws = round(fs/30)+1;% po = 21;% try%     yfil = sgolayfilt(yData, po, ws);% catch%     yfil = sgolayfilt(yData, po, ws+1); end hold on;plot( xData, yfil, '--g');hold off;
        xlim([-0.005, 0.020]);
        % ylim([-50 700]);
        drawnow;
        frame = getframe(gcf);
        writeVideo(videoObj, frame);
        % pause(0.5);
    end
    
    if 0
        figure(20);clf; 
        h = plot( xData-xData(1), yData );
        hold on;
        obj = feval(FitRes{1}, (xData));
        plot( xData-xData(1), obj, 'r');
        legend( h, 'Pressure vs. Time', 'Fitted model', 'Location', 'NorthEast', 'Interpreter', 'none' );
        xlabel( 'Time [s]', 'Interpreter', 'none' );
        ylabel( 'Foot pressure [Pa]', 'Interpreter', 'none' );
        grid on;
        
        title(sprintf([files(BIGIDX).name, ', trans num: ', num2str(jdx)]));
        
        ws = round(fs/30)+1;
        po = 21;
        try
            yfil = sgolayfilt(yData, po, ws);
        catch
            yfil = sgolayfilt(yData, po, ws+1);
        end
        plot( xData-xData(1), yfil, '--g');

        hold off;
        xlim([0.000, 0.250]);

        drawnow;
        pause(); 
    end
    
    if 0 
    fprintf('Pfoot_targ %1.1f, Beta %1.1f, Nu %1.5f\n',FitRes{jdx}.a,FitRes{jdx}.b,FitRes{jdx}.d);
    end

    % Store fit data
    % % % % % % % % % % % % % % % % % % % % % % % % % % afit(jdx)  = FitRes{jdx}.a; % Ptarg
    bfit(jdx)  = FitRes{jdx}.b; % beta
    dfit(jdx)  = FitRes{jdx}.d; % nu
    gofr2(jdx) = gof{jdx}.rsquare;

    
end % CLOSE LOOP over all transients of current file
fprintf(" (Done)\n");

Ptargfit = afit;
betafit  = bfit;
nufit    = dfit;

% Reduced Jet Velocity in S-S if Wm is known ==========================

RJV = sqrt(2*Pfoot_targ(:)/1.2)./(thedata.PR_params.Wm*f1(:));
Wm  = thedata.PR_params.Wm;

% %%%%%% SAVE processed data and results %%%%%%%%%%

datafilename = filename(1:end-4);

if 0
    save(['./processed/' datafilename '_PROCESSED.mat'],'f1','betafit','nufit',...
        'Area1','Area2','RJV','Wm','fs','PpalletB_targ','Pgroove_targ','t20groove',...
        'PRTgroove','Pfoot_targ','t20foot','PRTfoot','Ppipe_targ','t20mouth','PRTpipe',...
        'KeyMovingTime','gofr2','A2max_over_A1simult','A2max_over_A1target');
    
    fprintf('Ellapsed time with this file: %1.2f \n', toc(tinit));
end


end % BIGIDX of all files opened
close(videoObj);