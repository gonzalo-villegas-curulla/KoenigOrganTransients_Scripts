%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MODEL 03 pipe organ flow computation
%
% This script assumes a pressostat of p2Const [Pa]
% in the pallet-box. It computes the flow velocity from
% the pallet-box to the groove through the pallet-valve 
% with linear-ramp opening. It then computes the pressure 
% in the groove, the flow velocity from the groove to 
% the pipe foot inlet, the pressure in the foot, and the
% flow velocity out at the flue exit.

% The implementation consists of 3 unsteady Bernoulli's
% at the stages palletBox-groove, groove-foot, foot-mouth;
% in the cavities (groove and foot), mass conservation
% is applied under adiabatic assumption.

% The variables being solved in the ODE system are:
% y(1) = u3 = flow vel. from pallet box to valve-in-groove
% y(2) = p3 = pressure in the groove
% y(3) = u4 = flow vel. from groove to foot inlet
% y(4) = p4 = foot pressure
% y(5) = u5 = flue exit velocity
%        p5 = pressure @mouth (<p5> = 0)

% Gonzalo Villegas Curulla, Paris feb 2025.
% Running MATLAB R2021a. try() statements use the 
% wavelet toolbox and system identification toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear; fprintf('Starting script...\n');
try 
cd '/media/organ/ExtremeSSD/OrganPipe2023-2024/DataTransients/processed/NumModel'
end


% ========= Physical constants =======
rho = 1.2;
co = 340;
co2 = co^2;
P0 = 820;



% Vgrv = sth;
% Vf = sth_else;
% 
% A = Seff_pall *co2/Vgrv *sqrt(rho/P0);
% B = Seff_in   *co2/Vgrv *sqrt(rho/P0);
% C = Seff_in   *co2/Vgrv *sqrt(rho/P0);
% D = Seff_j    *co2/Vgrv *sqrt(rho/P0);
% 
% tau_pall = l_pall*sqrt(rho/P0);
% tau_in   = l_in  *sqrt(rho/P0);
% tau_j    = l_j   *sqrt(rho/P0);


PIPENUM = 10; TRANSNUM  = 20; 

files = dir('../../A*.mat');
load(fullfile('../../',files(PIPENUM).name));

x = PR_data.pressureData;
for idx = 1 : 6
    x(idx,:) = x(idx,:) - mean(x(idx,1:1e3));
end

% x(1,:); (Keller)
PRES  = x(2,:); 
PGRV  = x(3,:);
PF    = x(4,:);
PRAD  = x(5,:);
KEYV  = x(6,:);


[Lp,Vf,Pw,Tnhd,Sin,Sj,Hm,h,Wm,Rp,palletLHS,palletWid,palletRHS,palletHtraj] = getgeometry(PIPENUM);
Rin     = sqrt(Sin/pi);
Vgroove = Pw*0.51*0.05; % [m^3]
Sslot   = Pw*0.1298;

[KeyDownIdx,KeyUpIdx,KeyMovingTime,DurNotesInS, VelPeakIdxPos,VelPeakIdxNeg] = DetectVelocityPeaks_func(...
    KEYV, tvec, files(PIPENUM).name);

fs   = 51.2e3; 
dt   = 1/fs;
Tend = 0.150;
tmeas   = [0:dt:Tend]';
N    = numel(tvec);

PRECUT      = fix(0.010*fs);
SimulLength = fix(Tend*fs);
ll          = -PRECUT:-PRECUT + SimulLength ; % Mask, selection

pres = PRES(   KeyDownIdx(TRANSNUM)      + ll );
pgrv = PGRV(   KeyDownIdx(TRANSNUM)      + ll );
pf   = PF(     KeyDownIdx(TRANSNUM)      + ll );
prad = PRAD(   KeyDownIdx(TRANSNUM)      + ll );
keyv = abs(    KEYV(KeyDownIdx(TRANSNUM) + ll) ) ;
keyx = cumtrapz(dt, keyv);
keyx(PRECUT + (KeyUpIdx(TRANSNUM)-KeyDownIdx(TRANSNUM)) +fix(1e-3*fs): end) =  keyx(PRECUT + ( KeyUpIdx(TRANSNUM)-KeyDownIdx(TRANSNUM))+fix(1e-3*fs) -1) ;


% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%     DATA prepared: pres, pgrv, pf, prad, keyx
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


% ===================================
% Define Vena-Contracta 
% ===================================


VC3  = 0.6;  % pallet
VC4  = 0.42; % int
VC5  = 1.0;  % jet

ramptype = 'pall'; % 'rect'=linear, 'sin'=sinus, 'pall'=palletvalve

% ===================================
%   Physical constants and params
% ===================================

rho = 1.2;           % [kg/m^3]
c   = 340; c2 = c^2; % [m/s]

% ===================================
%   Organ sections
% ===================================


%%%%%%%%%%%%%%% (2)(Pallet box)  %%%%%%%%%%%%%%%

V2      = 0.186; % [m^3] Vres = Vbellow + Vtrunk + VpalletBox
p2Const = P0;    % [Pa] Pressostat value @ palletBox/Equiv.Resrv

%%%%%%%%%%%%%%% (3)(Groove)  %%%%%%%%%%%%%%%

GrvWidth  = Pw; %10.9e-3; % (6.2-15.9)mm
GrvLen    = 0.510 ;
GrvHeight = 0.0496; 

V3        = 1.3*GrvWidth*GrvLen*GrvHeight;   % [m^3] Groove volume 
l3        = 1.0*8e-3;     % [m] PalletValve-Groove slot inlet

% S3 "slot" measured: L=129.8mm, Width=[6.2,15.9]mm
Spallet   = 0.1298*GrvWidth; % [m^2]


% +++++  Create a valve opening ramp:  ++++++++++++
ValveOpenInit = 0.100; % [s] T-Start opening-ramp pallet valve
ValveOpenEnd  = 0.130; % [s] T-Finish opening-ramp (DROPIC robot time)


S3                    = zeros(N,1);
S3(tvec>ValveOpenEnd) = Spallet;
switch ramptype
    case 'rect'
        ll     = find( ValveOpenInit<=tvec & tvec<=ValveOpenEnd );
        ramp   = linspace(0,Spallet,numel(ll));
        S3(ll) = ramp;
    case 'sin'
        ll     = find( ValveOpenInit<=tvec & tvec<=ValveOpenEnd );
        sintime = linspace(-pi/2,pi/2, length(ll));
        vec = 0.5+ 0.5*sin(sintime);
        vec = Spallet*vec;
        S3(ll) = vec;
    case 'pall'
        % Spallet*[0, ... ,1] 
        
        % (1):
        % S3 = (Pw + 0.1298) * keyx;
        
        % (2):
        faco = 1.5 + 0*0.999;
        S3 = (faco*min(palletLHS,palletHtraj)+ palletWid+faco*min(palletRHS,palletHtraj))*keyx/max(keyx);
        
        faco = 0.5;
        S3 = (faco*min(palletLHS,palletHtraj)+ palletWid+faco*min(palletRHS,palletHtraj))*keyx/max(keyx);
        
        % (3):
        % S3 = keyx*Sslot;
end 
mav    = dsp.MovingAverage(fix(0.001*fs));
% S3  = mav(S3); % Optionally, smooth the ramp



%%%%%%%%%%%%%%% (4)(Pipe foot)  %%%%%%%%%%%%%%%

V4  = Vf;      
l4  = 0 + 2e-3;  
Sin = pi*Rin^2;  


%%%%%%%%%%%%%%% (5)(Flue exit)  %%%%%%%%%%%%%%% 

l5  = 0 + 2e-3; 
S5  = Sj ;


% ===================================
%  APPLY VENA CONTRACTA where relevant
% ===================================

S3vec = S3      *VC3;
S4    = Sin     *VC4;
S5    = Sj      *VC5; 


% ===================================
%     Pre-compute coefficients
% ===================================

params = [rho l3 c2 V3 S4 l4 V4 S5 l5, p2Const];


% ===================================
% Solve system 
% ===================================
fprintf('Starting ODE solver...\n');
tic

IC     = 1e-6*[1,1,1,1,1]; % Avoid null, solver takes longer.
tstart = 0;
tfinal = Tend;
y0     = IC;
refine = 4;
tcount = 0;
ctr    = 0;

try 
    opts  = odeset('RelTol',1e-5,'AbsTol',1e-5,'Events',@events1,'Events',@events2,'Refine',refine,'OutputFcn',@odetpbar);
catch 
    opts  = odeset('RelTol',1e-5,'AbsTol',1e-5,'Events',@events1,'Events',@events2,'Refine',refine);
    %%% opts  = odeset('RelTol',1e-5,'AbsTol',1e-5,'Events',@events1,'Refine',refine,'OutputFcn',@odetpbar);
end


tsimul  = tstart;
yout  = y0;
teout = [];
yeout = [];
ieout = [];

while tcount < Tend

% Solvers: 15s, 113 (accurate, sometimes slow), 23, 23s
[t_ode,y,te,ye,ie] = ode113(@(t_ode,y) myodes(t_ode, y , tmeas, S3vec, pres,  params), [tstart tfinal], y0, opts); 

   % Accumulate output. If solution goes to zero, re-start with its
   % perturbated conditions.
   nt    = length(t_ode);
   tsimul  = [tsimul; t_ode(2:nt)];
   yout  = [yout; y(2:nt,:)];
   teout = [teout; te]; 
   yeout = [yeout; ye];
   ieout = [ieout; ie];
   
   try
   y0(1) = ye(1);
   y0(2) = ye(2);
   y0(3) = ye(3);
   y0(4) = abs( ye(4) )+0.5;
   y0(5) = abs( ye(5) )+0.001;
   end
   
   % Guess of a valid first timestep: length of the last valid timestep
   % 'refine' defaulted to 4
   opts = odeset(opts,'InitialStep',t_ode(nt)-t_ode(nt-refine),'MaxStep',t_ode(nt)-t_ode(1));

   tstart = t_ode(nt);
   tcount = t_ode(end);
   ctr    = ctr+1;
end

% ===================================
%       Read-out solutions
% ===================================
fprintf('(done)\n');
toc


% ===================================
%       Fit logistic in Pfoot
% ===================================

% PREPARE FIT OPTIONS ===========================
% Set up fittype and options.
func = fittype( 'a./(1+d*exp(-b*(x-c))).^(1/d)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display    = 'Off';
opts.Robust     = 'Bisquare'; % LAR, Off, Bisquare

guessP = mean(yout(end-100:end,4)); % Last values of foot pressure (~SS)

% %                Ptarg  Beta  t_o   Nu
opts.Lower      = [1      1     1e-3  1e-3];
opts.Upper      = [1000   15e3  3     30   ];
opts.StartPoint = [guessP 100   0.050 1];
    
% xData = tsimul; % Time stamps after solver
% yData = yout(:,4); % Solved foot pressure
xData = tmeas(:);
yData = pf(:); % Solved foot pressure
[FitRes, gof] = fit( xData, yData, func, opts );
    
afit  = FitRes.a;
bfit  = FitRes.b;
cfit  = FitRes.c;
dfit  = FitRes.d;
gofr2 = gof.rsquare;

% Fitted func result:
% fittedPfoot = afit./(1+exp(-bfit*(tsimul-cfit))).^(1/dfit);
fittedPfoot = afit./(1+exp(-bfit*(tmeas-cfit))).^(1/dfit);


% ===================================
% Resample simulated results to measured query points
% ===================================
u3 = yout(:,1);
p3 = yout(:,2);
u4 = yout(:,3);
p4 = yout(:,4);
u5 = yout(:,5);

p3 = interp1(tsimul, p3, tmeas);
p4 = interp1(tsimul, p4, tmeas);

% ===================================
% Find and shift key->foot delay
% ===================================
ThR = 0.5*p3(end);
shif = find(p3>ThR,1,'first') - find(pgrv>ThR,1,'first');
p3 = circshift(p3, -shif);
p3(1:abs(shif)) = pf(abs(shif)+1);

ThR = 0.5*p4(end);
shif = find(p4>ThR,1,'first') - find(pf>ThR,1,'first');
p4 = circshift(p4, -shif);
p4(1:abs(shif)) = pf(abs(shif)+1);

% ===================================
%           Plot results
% ===================================




FSZ = 14;
figure(1); clf;

ax(1) = subplot(411);
plot(tmeas, pres);
grid on; box on; ylabel('P PalletBox','fontsize',FSZ);
yyaxis right;
plot(tmeas, S3vec,'r');ylabel('Valve area','fontsize',FSZ);

ax(2) = subplot(412);
plot(tmeas,   pgrv); hold on;
plot(tmeas, p3 );
grid on; box on; ylabel('P Groove','fontsize',FSZ);

ax(3) = subplot(413);
plot(tmeas,   pf); hold on;
plot(tmeas, p4);
% plot(tmeas, fittedPfoot,'-m');
grid on; box on; ylabel('P Foot','fontsize',FSZ);
xlabel('Time ');

ax(4) = subplot(414);
plot(tmeas, prad);
% hold on;
grid on; box on; ylabel('P Rad','fontsize',FSZ);

% ax(4) = subplot(514);
% plot(tsimul, p4, '-o');
% hold on;
% plot(tsimul, fittedPfoot);
% grid on; box on; ylabel('PRESS. Foot','fontsize',FSZ);
% title(sprintf('Beta= %1.2f; Nu: %1.3f', bfit,dfit));
% %
% ax(5) = subplot(515);
% plot(tsimul, u5, '-o');
% grid on; box on; ylabel('JET VELOC.','fontsize',FSZ); xlabel('time [s]','fontsize',FSZ);


linkaxes(ax,'x');
xlim([0 Tend]);


%%
FSZ = 14;
figure(1); clf;
ax(1) = subplot(511);
plot(tsimul, u3, '-o');
grid on; box on; ylabel('FLOW Pallet','fontsize',FSZ);
yyaxis right;
plot(tmeas, S3vec,'r');ylabel('Valve area','fontsize',FSZ);
%
ax(2) = subplot(512);
plot(tsimul, p3, '-o');
grid on; box on; ylabel('PRESS. Groove','fontsize',FSZ);
%
ax(3) = subplot(513);
plot(tsimul, u4, '-o');
grid on; box on; ylabel('FLOW Foot IN','fontsize',FSZ);
%
ax(4) = subplot(514);
plot(tsimul, p4, '-o');
hold on;
plot(tsimul, fittedPfoot);
grid on; box on; ylabel('PRESS. Foot','fontsize',FSZ);
title(sprintf('Beta= %1.2f; Nu: %1.3f', bfit,dfit));
%
ax(5) = subplot(515);
plot(tsimul, u5, '-o');
grid on; box on; ylabel('JET VELOC.','fontsize',FSZ); xlabel('time [s]','fontsize',FSZ);


linkaxes(ax,'x');
xlim([0 Tend]);


figure(2); clf; 
plot(tmeas, pf);
hold on;
plot(t_ode, p4);

%%

% Below here, not very relevant calculations:
try
    % Spectrum of a signal (freqs. due to geometry changes):
    F = griddedInterpolant(tsimul,p4); %To choose from: u3,p3,u4,p4,u5
    X = F(tvec);
    figure(2); clf;
    Wlen = 500;
    OL = fix(Wlen*0.75);
    Nfft = 2^15;
    cwt(X,[],fs)

    % Attemp modal freq extraction from the pallet valve ramp function:
    win = hann(length(X));
    [frf,f] = modalfrf(S3vec',X',fs,win,'Sensor','vel'); %acc,vel,dis
    [ModalNatFreq,DampRatio] = modalfit(frf,f,fs,1,'FitMethod','lsce'); %lsce,lsrf,pp
    if 0
            figure(3); clf;
            modalfit(frf,f,fs,1,'FitMethod','lsce')   
    end
end

% ==================================
% %         FUNCTIONS 
% ==================================


function dydt = myodes(t,y,tvec,S3vec,pres, params)

rho = params(1);
l3  = params(2);
c2  = params(3);
V3  = params(4);
S4  = params(5);
l4  = params(6);
V4  = params(7);
S5  = params(8);
l5  = params(9);
p2Const = params(10);

S3interp = interp1(tvec, S3vec, t, 'next');% nearest, next, pchip, spline, 
p2int    = interp1(tvec, pres, t, 'next');


dydt = zeros(5,1);
P5   = 0;

    if (S3interp~=0)
%         dydt(1) = 1/rho/l3*(p2Const  -y(2))  -1/2/l3*y(1)^2    ;
        dydt(1) = 1/rho/l3*(p2int  -y(2))  -1/2/l3*y(1)^2    ;
        
        dydt(2) = rho*c2/V3* ( y(1)*S3interp  -S4*y(3) )  ;
        dydt(3) = 1/rho/l4*(y(2)  - y(4))  -1/2/l4*y(3)^2 ;
        dydt(4) = rho*c2/V4*( y(3)*S4  -S5*y(5) )         ;
        dydt(5) = 1/rho/l5*(y(4)  - P5)  -1/2/l5*y(5)^2   ;
    end

    
    

    % dydt(1) = 
    % dydt(2) = 
    % dydt(3) = 
    % dydt(4) = 
    % dydt(5) = 


end

% Zero detection of solution event triggers:
function [value,isterminal,direction] = events1(t_ode,y)
value = y(4);     % detect height = 0
isterminal = 1;   % stop the integration
direction = 0;   % negative direction
end

function [value,isterminal,direction] = events2(t_ode,y)
value = y(5) | y(4);     % detect height = 0
isterminal = 1;   % stop the integration
direction = 0;   % negative direction
end
