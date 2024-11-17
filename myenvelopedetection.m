noisetime = 0.100;
NumPers   = 15; % Number periods for target averaging
T1estim   = 1 / f1estim; 
hL   = fix(0.5*length(pipecurr)); % half length of the signal
L    = 2*hL;
time = (0:hL-1) / fs;
fact = 0.9;

%%%%%%%%%%%%%%%%%%%%%%%%%
% TOTAL PRESSURE 
%%%%%%%%%%%%%%%%%%%%%%%%%

signal = pipecurr(1:hL);     % Keep 1st half
signal(time<=noisetime) = 0; % delete noise initial

% At the end of half-signal, find target amplitude using NumPers periods
target_amp = mean(maxk(signal(end-fix(fs*NumPers*T1estim/1):end),NumPers)); 
TOTAL_target_amp = target_amp;
separ      = fix(0.9*fs*T1estim);


% Ballpark beginning of precursor (get its MaxNeg from the Fundamental)
tmp = COMP(1,1: fix(   fs*(noisetime+5*T1estim) ) );
[~,maxnegprecursor_idx] = max(-tmp);
% Let's say precursor start at the last positive sample before negative max
init_precursor_idx = find(signal(1:maxnegprecursor_idx)>0,1, 'last' );
if isempty(init_precursor_idx)
    init_precursor_idx = fix(fs*noisetime);
    % Otherwise, run gvcstft and find where energy shoots for the precursor
end


% Create masks for (init noise)(prominentTransient until 90%)(preSteady&Steady)
mask1 = 1:init_precursor_idx;
mask2 = mask1(end)+1: find( max(0,signal)/target_amp>0.9, 1, 'first' ); 
mask3 = mask2(end)+1:hL;

% /!\ Noise init 
signal(mask1) = 0; % Delete noise until beginning of precursor

% /!\ growing
wid = fix( (0.6)*(1/2)*(1/3)*T1estim*fs); % (Headroom)(Half-period=peak)(shortest Period I careabout)
[pks2,locs2]  = findpeaks(signal(mask2), 'MinPeakDistance',separ, 'MinPeakWidth',wid  );
    pks2  = [0,pks2];locs2 = [1,locs2]; % add the beginning of the precursor for interpolation (avoid initial spike)
envel2 = interp1(time(mask1(end)+locs2),pks2, time(mask2), 'pchip');


% /!\ steady-ish
[pks3,locs3]    = findpeaks(signal(mask3),'MinPeakHeight',0.6*target_amp,'MinPeakDistance',separ);
envel3          = interp1(  time(mask2(end)+locs3), pks3, time(mask3), 'pchip');

% Recompose whole envelope in (0,hL)
enveltot        = zeros(size(signal));
enveltot(mask2) = envel2;
enveltot(mask3) = envel3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNDAMENTAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

signal      = COMP(1,1:hL);
signal(time<=noisetime) = 0;
envel_first = zeros(size(signal));

% separ = fix(0.9*fs*T1estim/1);
target_amp = mean(maxk(signal(end-fix(fs*NumPers*T1estim/1):end),NumPers)); 

% mask1 = same (precursor)
mask2 = mask1(end)+1: find( max(0,signal)/target_amp>0.9, 1, 'first' ); 
mask3 = mask2(end)+1:hL;

signal(mask1) = 0;
wid           = fix( (0.6)*(1/2)*(1/3)*T1estim*fs); 
[pks2,locs2]  = findpeaks(signal(mask2), 'MinPeakDistance',separ, 'MinPeakWidth',wid  ); pks2  = [0,pks2];locs2 = [1,locs2]; 
envel2        = interp1(time(mask1(end)+locs2),pks2, time(mask2), 'pchip');
[pks3,locs3]  = findpeaks(signal(mask3),'MinPeakHeight',0.6*target_amp,'MinPeakDistance',separ);
envel3        = interp1(  time(mask2(end)+locs3), pks3, time(mask3), 'pchip');

envel_first(mask2) = envel2;
envel_first(mask3) = envel3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECOND HARMONIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

signal       = COMP(2,1:hL);
signal(time<=noisetime) = 0;
envel_second = zeros(size(signal));

separ = fix(0.95*fs*T1estim/1); % it should be T1/2, but this causes bleeding
target_amp = mean(maxk(signal(end-fix(fs*NumPers*T1estim/2):end),NumPers)); 
target_second = target_amp;

% mask1 = same 
mask2 = mask1(end)+1: find( max(0,signal)/target_amp>0.9, 1, 'first' ); 
mask3 = mask2(end)+1:hL;

signal(mask1) = 0;
wid           = fix( (0.6)*(1/2)*(1/3)*T1estim*fs); 
[pks2,locs2]  = findpeaks( max(0,signal(mask2)), 'MinPeakDistance',separ, 'MinPeakWidth',wid  ); pks2  = [0,pks2];locs2 = [1,locs2]; 
envel2        = interp1(time(mask1(end)+locs2),pks2, time(mask2), 'pchip');
[pks3,locs3]  = findpeaks(signal(mask3),'MinPeakHeight',0.6*target_amp,'MinPeakDistance',separ);
envel3        = interp1(  time(mask2(end)+locs3), pks3, time(mask3), 'pchip');

envel_second(mask2) = envel2;
envel_second(mask3) = envel3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIRD HARMONIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

signal       = COMP(3,1:hL);
signal(time<=noisetime) = 0;
envel_third = zeros(size(signal));

separ = fix(0.9*fs*T1estim/1); % This should be T1/3, but it causes bleeding
target_amp = mean(maxk(signal(end-fix(fs*NumPers*T1estim/3):end),NumPers)); 

% mask1 = same 
mask2 = mask1(end)+1: find( max(0,signal)/target_amp>0.9, 1, 'first' ); 
mask3 = mask2(end)+1:hL;

signal(mask1) = 0;
wid           = fix( (0.6)*(1/2)*(1/3)*T1estim*fs); 
[pks2,locs2]  = findpeaks(signal(mask2), 'MinPeakDistance',separ, 'MinPeakWidth',wid  ); pks2  = [0,pks2];locs2 = [1,locs2]; 
envel2        = interp1(time(mask1(end)+locs2),pks2, time(mask2), 'pchip');
[pks3,locs3]  = findpeaks(signal(mask3),'MinPeakHeight',0.6*target_amp,'MinPeakDistance',separ);
envel3        = interp1(  time(mask2(end)+locs3), pks3, time(mask3), 'pchip');

envel_third(mask2) = envel2;
envel_third(mask3) = envel3;

envel_sum = envel_first+envel_second+envel_third;