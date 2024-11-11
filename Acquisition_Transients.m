%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Acquisition Koenig jussieu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc;

%----------------------------------------------------
%          0. Intro
%----------------------------------------------------

 % -----------------------------
PR_metadata.filenamePrepend = 'A07_pitchFsharp_';%
PR_metadata.extraInfo       = '';
PR_metadata.savePath        = ['']; % Existing folder in path
PR_params.Tacq              = 90+0*300; % Acquisition time [s]
% -----------------------------
PR_params.SR         = 51200;     % Sampling frequency requested [Hz]
PR_params.chanNumber = 6;         % Number of channels
% -----------------------------
PR_metadata.Temperature = 20.3;
PR_metadata.Humidity    = 43.3;
% -----------------------------

PR_metadata.mics = 'Keller_11337_11545_11312_3186627_Vibro';

% Probe sentitivities 
PR_params.sens(1)     = 1/0.0047; % Keller
PR_params.sens(2)     = 1/0.001; % Endevco 11337
PR_params.sens(3)     = 1/0.001; % Endevco 11545
PR_params.sens(4)     = 1/0.001; % Endevco
PR_params.sens(5)     = 1/0.00316; % BK
PR_params.sens(6)     = 0.125; % % Vibro



%----------------------------------------------------
%          1. Setup session
%----------------------------------------------------

daqreset
PR_metadata.device         = daq.getDevices;
Session                    = daq.createSession('ni');
Session.Rate               = PR_params.SR;
Session.DurationInSeconds  = PR_params.Tacq;

ch = addAnalogInputChannel(Session, PR_metadata.device(1).ID, [0:3] ,'Voltage');
% ch = addAnalogInputChannel(Session, PR_metadata.device(2).ID, [0:1] ,'Voltage'); 



  for idx = 1 : length(ch)
     ch(idx).Coupling = 'DC';%'DC';
  end
 SR = Session.Rate;
 
%----------------------------------------------------
%          2. Acquisition
%----------------------------------------------------
fprintf('Acquisition starting...\n');
[data,tvec, what] = startForeground(Session);
fprintf('Acquisition finished!\n');

%Correct amplitude sensitivities of each transducer:
for idx = 1 : numel(PR_params.sens)
    PR_data.pressureData(idx,:) = data(:,idx)*PR_params.sens(idx);
end
PR_data.Time         = tvec;
PR_params.SRadjusted = Session.Rate;
  
%----------------------------------------------------
%          3. options.plots
%----------------------------------------------------

LW = 1;
figure(1); clf;
hold on;
for idx = 1 : (min(size(PR_data.pressureData))-1)
   plot(tvec, (PR_data.pressureData(idx,:)-mean(PR_data.pressureData(idx,1:100))), 'linewidth',LW);
   hold on;
end
grid on;
xlabel('time [s]','fontsize',18);
ylabel('Pressure','fontsize',18);
yyaxis right;
plot(tvec, (PR_data.pressureData(end,:)-mean(PR_data.pressureData(end,1:100))), 'linewidth',LW);
ylabel('key velocity','fontsize',18);
legend('bellow','palletB','groove','foot','mouthRad','keyVeloc.','fontsize',18)

% % %% FOR printing the peak2peak/2 RMS [ =/sqrt(2)] SPL when calibrating with B&K's Pistonphone
% % fprintf('Peak to Peak: %1.3f dB\n',db(peak2peak(PR_data.pressureData)/2/sqrt(2)/20e-6));


% %======================================
% pfoot = PR_data.pressureData(4,:);
% ll    = fix(length(pfoot)/2): fix(length(pfoot)/2) + fix(SR*0.500);
% avpfoot = mean(pfoot(ll));
% avuj  = sqrt(2*avpfoot/1.2);
% Wm    = 7.81e-3;
% pmouth = PR_data.pressureData(5,:);
% R = yin(pmouth(ll)',SR);
% f1estim = mean(R.f0_Hz,'omitnan');
% 
% mytheta = avuj/(Wm*f1estim);
% fprintf('Calc RJV: %1.3f \n',mytheta)
% %======================================


%----------------------------------------------------
%          4. Save data
%----------------------------------------------------

tempInfoDevice = PR_metadata.device;
PR_metadata    = rmfield(PR_metadata, 'device');
PR_metadata.device.Model       = tempInfoDevice.Model;
PR_metadata.device.ID          = tempInfoDevice.ID;
PR_metadata.device.Description = tempInfoDevice.Description;

thisDate = datetime('now', 'format', 'yyyyMMdd_HHmmss');
thisDate = char(thisDate);
PR_metadata.dataFilename = [PR_metadata.filenamePrepend ...
                thisDate '.mat'];

save([pwd filesep PR_metadata.savePath filesep PR_metadata.dataFilename],...
    '-v7.3', '-nocompression');

% % % %Save fig
% % % figname = PR_metadata.dataFilename;
% % % figname(end-2:end) = 'fig';
% % % savefig(gcf, figname);

fprintf('Data saved in %s\n', PR_metadata.dataFilename);

% ===========================================

% run DetectVelocityPeaks.m