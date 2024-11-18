try
    cd /run/media/gvc/ExtremeSSD/OrganPipe2023-2024/DataTransients/
end 

clc, clear;
addpath('./processed/');
% ============= LOADS =========================

thepath = './processed/';
thepath = './processed/FIT_with_plus3PRT/';
filecard = 'A*PROCESSED.mat';
files = dir([thepath,filecard]);


% Resoantor length: (56)
load('./processed/Geometry/Lp_m.mat'); 
LP = Lp_m; clear Lp_m

% Foot volume: (56)
load('./processed/Geometry/Vf_m3.mat'); 
VF = Vf_m3; clear Vf_m3

% Width of pallet window slot (56)
load('./processed/Geometry/KoenigPalletWindWidth_m.mat');
PW    = palletwinwidth_m; clear palletwinwidth_m

% Tone hole diam, on the wooden top board: (56)
load('./processed/Geometry/ToneHoleDiam_m.mat');
TNHD = ToneHoleDiam; clear ToneHoleDiam

% Toe hole area (foot inlet): (56)
load('./processed/Geometry/InletManipFitted_mm2.mat'); 
INLET = InletManipFitted_mm2* 1e-6; clear InletManipFitted_mm2

% Flue exit area (S_jet): (56)
load('./processed/Geometry/FlueExitManipFitted_mm2.mat');  
SJET    = FlueExitManipFitted_mm2 * 1e-6; clear FlueExitManipFitted_mm2

% Flue exit little height: (56)
load('./processed/Geometry/SmallhManipFitted_mm.mat');
h     = hManipFitted_mm * 1e-3; clear hManipFitted_mm

%FluExit width, big Hm: (56)
load('./processed/Geometry/BigHManipFitted_mm.mat'); 
H     = HManipFitted_mm * 1e-3; clear HManipFitted_mm

% Mouth cutupt distance to labium (56):
load('./processed/Geometry/WmManipFitted_mm.mat');
WM    = WmManipFitted_mm * 1e-3; clear WmManipFitted_mm 

% Resonator diameter (56)
load('./processed/Geometry/DpManipFitted_mm.mat');
DP    = DpManipFitted_mm*1e-3; clear DpManipFitted_mm

% F1 (averaged) for the measured pipes (16) <<<<!
load('./processed/Geometry/f1ManipMean.mat');
% >> F1MEAN

load('./processed/MaxNegKeyVelocities.mat');
KeyVel = abs(MAXNEGVELOS);

% ============= PARAMETERS =========================

rho   = 1.2;
c     = 340;
maskpipes  = [3,4,5,6,7,9,10,11,13,15,17,19,24,25,27, 29,32, 34,37,39,41,44]';



PalletDepth = 1e-3*129.8; %  [m]
PALLAREA    = PalletDepth * PW;  % Slot area covered by the valve [m^2]
VGROOVE     = PW *0.51*0.05;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefs and params 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vena contractai:
VCpallet = 0.62;
VCinlet  = 0.62;
VCjet    = 0.95;
VCpallet=1;VCinlet=1;VCjet=1;

NumTransMax = 100;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Q-factors for fundamentals
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RADS     = 0.5*DP;
RADSmask = RADS(maskpipes);
DPmask   = DP(maskpipes);
LPmask   = LP(maskpipes);

n = 1; % mode
chi       = 3e-5 * sqrt(n*F1MEAN)./RADSmask;
OM        = n*pi*c./(LPmask + 1.2*RADSmask);
part1     = OM/(2*c);
part2     = (LPmask + 1.2*RADSmask)./(chi.*LPmask + (OM.^2.*RADSmask.^2)/(2*c^2) )  ;
QFAC1     = part1.*part2;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MX            = nan*ones(NumTransMax, length(files), 49);
vecmeanfreqs  = zeros(length(files),1);
PRTMX         = nan*ones(NumTransMax, length(files), 3);
PRTpipeMedian = zeros(length(files),1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Populate data in matrices 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for idx = 1 : length(files)
    
    filename = files(idx).name;
    try
        load(filename);
    catch
        load([thepath,filename]);
    end
    
    numpipe    = str2num(filename(2:3));    % MODIFY IF NECSSARY<<<<<<<<<<
    
    
    % #56 data:
    oneLP      = LP(numpipe); % Pipe length
    oneVF      = VF(numpipe); % Volume foot 
    onePW      = PW(numpipe); % Pallet valve slot width 
    oneTNHD    = TNHD(numpipe);
    oneINLET   = INLET(numpipe);
    oneSJET    = SJET(numpipe);
    oneh       = h(numpipe);
    oneH       = H(numpipe); 
    oneWM      = WM(numpipe); % Mouth W cutup
    oneDP      = DP(numpipe); % Diameter pipe
    oneVGROOVE = VGROOVE(numpipe);
    
    % #16 data:
    onef1      = F1MEAN(idx);
    oneQFAC1   = QFAC1(idx);
    
    %%%%%%%%%% GEOMETRICAL ASPECTS %%%%%%%%%% 
%     NumTransMax = 1; % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    MX(:, idx, 1)    = oneLP*ones(NumTransMax,1);
    MX(:, idx, 2)    = oneVF*ones(NumTransMax,1);
    MX(:, idx, 3)    = onePW*ones(NumTransMax,1);
    MX(:, idx, 4)    = oneTNHD*ones(NumTransMax,1);
    MX(:, idx, 5)    = oneINLET*ones(NumTransMax,1);
    MX(:, idx, 6)    = oneSJET*ones(NumTransMax,1);
    MX(:, idx, 7)    = oneh*ones(NumTransMax,1);
    MX(:, idx, 8)    = oneH*ones(NumTransMax,1);
    MX(:, idx, 9)    = oneWM*ones(NumTransMax,1);
    MX(:, idx, 10)   = oneDP*ones(NumTransMax,1);
    MX(:, idx, 11)   = oneVGROOVE*ones(NumTransMax,1);
    MX(:, idx, 12)   = oneQFAC1*ones(NumTransMax,1);
    
    
    
    %%%%%%%%%% S-S ASPECTS %%%%%%%%%% 
    
    % S1/S2 factor, Spall/SJET(Pf,Pgrv) Remplissage factor:
    Remplissage = real(sqrt(Pfoot_targ ./ (Pgroove_targ-Pfoot_targ ) ) ); 
    Qpall2grv   = VCpallet * PALLAREA(numpipe) * real(sqrt( (PpalletB_targ - Pgroove_targ) *2/rho));
    Qgrv2ft     = VCinlet  * INLET(numpipe)    * real(sqrt( (Pgroove_targ  - Pfoot_targ  ) *2/rho));
    Qjet        = VCjet    * SJET(numpipe)     * real(sqrt( (Pfoot_targ    - 0           ) *2/rho));
    
    MX(find(f1)                 , idx, 13) = f1(:);
    MX(find(RJV)                , idx, 14) = RJV;
    MX(find(Qpall2grv),           idx, 15) = Qpall2grv;
    MX(find(Qgrv2ft)  ,           idx, 16) = Qgrv2ft;
    MX(find(Qjet),                idx, 17) = Qjet;
    MX(find(Remplissage)       ,  idx, 18) = Remplissage(Remplissage~=0);
    
    
    MX(find(Pfoot_targ), idx, 19) = PpalletB_targ;
    MX(find(Pfoot_targ), idx, 20) = Pgroove_targ;
    MX(find(Pfoot_targ), idx, 21) = Pfoot_targ;
    MX(find(Pfoot_targ), idx, 22) = Ppipe_targ;
    
    
    if 0
        MX(find(Pgroove_targ./PpalletB_targ), idx, 23)  = Pgroove_targ./PpalletB_targ;
        MX(find(Pfoot_targ./Pgroove_targ)   , idx, 24)  = Pfoot_targ  ./Pgroove_targ;
        MX(find(Ppipe_targ./Pfoot_targ)     , idx, 25)  = Ppipe_targ  ./Pfoot_targ;
    else
        MX(find(Pgroove_targ), idx, 23)     = Pgroove_targ - PpalletB_targ;
        MX(find(Pfoot_targ), idx, 24)       = Pfoot_targ   - Pgroove_targ;
        MX(find(Ppipe_targ), idx, 25)       = Ppipe_targ   - Pfoot_targ;
    end
    
    %%%%%%%%%% TRANSIENT ASPECTS %%%%%%%%%% 
    
    t20foot=t20foot(:);t20groove=t20groove(:);t20mouth=t20mouth(:);
    PRTgroove = PRTgroove(:);PRTfoot=PRTfoot(:);PRTpipe=PRTpipe(:);
    
    MX(find(betafit)  , idx, 26)  = betafit(:); % Foot pressure transient slope
    MX(find(nufit)    , idx, 27)  = nufit(:);   % Foot pressure transient onset sharpness
    
    
    MX(find(PRTgroove), idx,28) = PRTgroove;
    MX(find(PRTfoot)  ,idx,29)  = PRTfoot;
    MX(find(PRTpipe)  ,idx,30)  = PRTpipe;
    
    MX(find(PRTfoot)     ,idx, 31) = PRTfoot./PRTgroove;
    MX(find(PRTpipe)     ,idx, 32) = PRTpipe./PRTfoot;
    
    MX(find(t20groove), idx, 33) = t20groove(t20groove~=0);
    MX(find(t20foot),   idx, 34) = t20foot(t20foot~=0);
    MX(find(t20mouth),  idx, 35) = t20mouth(t20mouth~=0);
    
    MX(find(t20foot-t20groove), idx, 36) = (t20foot-t20groove);
    MX(find(t20mouth-t20foot),  idx, 37) = (t20mouth-t20foot);
    
    MX(find(betafit)            , idx, 38)   = Area1(:);
    MX(find(betafit)            , idx, 39)   = Area2(:);
    
    
    MX(:           , idx, 40)   = ones(NumTransMax,1)*(INLET(numpipe)/PALLAREA(numpipe));
    MX(:           , idx, 41)   = ones(NumTransMax,1)*(SJET(numpipe)/INLET(numpipe));
    
    MX( : , idx, 42)            = KeyVel(:,idx);

    MX(find(A2max_over_A1simult),idx,43)       = A2max_over_A1simult;
    MX(find(A2max_over_A1target),idx,44)       = A2max_over_A1target;
    MX(find(A2max_over_A2target),idx,45)       = A2max_over_A2target;
    MX(find(pf_at_a2max - pm_at_a2max),idx,46) = pf_at_a2max - pm_at_a2max; % DeltaP (foot-mouth)
    MX( find(a2max_vec),idx,47)                = a2max_vec;
    MX(find(max_a2_over_a1), idx, 48)          = max_a2_over_a1;
    MX(find(gofr2), idx, 49)                   = gofr2;
    
    
    
   
    %%%%%%%%%% MARGINAL ASPECTS %%%%%%%%%% 
    PRTpipeMedian(idx) = median(PRTpipe(PRTpipe~=0),'omitnan');
    vecmeanfreqs(idx)  = mean(f1);
    freqref            = 440;
    freqlogax          = 12*log2(vecmeanfreqs/freqref);
    
    
end    
%%%%%%%%%%%%%%%%%%%% PREPARE PCA()  %%%%%%%%%%%%%%%%%%%%



MXresh = nan*ones(size(MX,1)*size(MX,2),size(MX,3));
for idx = 1 : size(MX,3)
   MXresh(:,idx) = reshape( MX(:,:,idx),numel(MX(:,:,idx)),1 ) ;  
end
    


% [1]:Lp        [2]:Vf       [3]:PWidth  [4]:ToneHoleDiam  [5]:Inlet 
% [6]:Sjet      [7]:h        [8]:H       [9]:Wm            [10]:Dpipe        
% [11]:Vgroove  [12]:Qfactor        

% [13]:f1          [14]:theta        [15]:Qpall2groove   
% [16]:Qgrv2foot   [17]:Qjet         [18]:Remplissage  
% [19]:Ppall targ  [20]:Pgrove targ  [21]:Pfoot targ  [22]:Prad targ

% [23]:Ptarg grv-pall       [24]:Ptarg foot-grv       [25]:Ptarg rad-foot

% [26]:beta                 [27]:nu   

% [28]:PRTgrv               [29]:PRTfoot              [30]:PRTrad
% [31]:PRTfoot/grv          [32]:PRTrad/foot  
% [33]:t20 groove           [34]:t20 foot             [35]:t20 rad
% [36]:t20 foot-groove      [37]:t20 rad-foot
% [38]:Area1                [39]:Area2
% [40]:Sin/Spall            [41]:Sjet/Sin             [42]: KeyVel

% [43]: A2max_over_A1simult [44]: A2max_over_A1target [45]: a2max_over_a2target
% [46]: DeltaP(foot-mouth)_at_a2max
% [47]: a2max_vec        [48]: max_a2_over_a1 (after smooth)\in(t20_f,t80_f+50PRT)
% [49]: gofr2 (r-squared goodness of logistic fit)



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Steady-State analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho = 1.2;

      %   S_in                     SS target    SS target
      %   geom                      Preserv      Pfoot 
Qin = 1*MX(:,:,5) .* sqrt( 2/rho*( MX(:,:,19) - MX(:,:,21) ) );

       % Geom                      SS target   SS target
       % S_j                        P_foot      P_mouth
Qj  = 1*MX(:,:,6).*sqrt( 2/rho * ( MX(:,:,21) -      0     ) );

figure();
                % All the f1's
scatter( 12*log2(  MX(:,:,13   )/440) , abs(Qin-Qj)./Qin ,'b','filled');
% scatter( 12*log2(  MX(:,:,13   )/440) , Qj./Qin ,'b','filled');
grid on; xlabel('tessitura');ylabel('Qj/Qin'); box on;
% ylim([0 1]);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Transient analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assess the goodness of fit for the cases: t80+1PRT, t80+2PRT, t80+3PRT
figure();
scatter( 12*log2(MX(:,:,13)/440), MX(:,:,49), 'b', 'filled');
xlabel('tessitura semitones');
ylabel('R-squared goodness of fit');
grid on; box on;



%% max(a2/a1) during transient
% butter(4), filtfilt() of envel_first and envel_second
% then search max within t20_f and t80_f+50*PRT

mask = [1:22]; %mask(7)=[]; % >>> If you want to remove a single dataset

figure();
scatter( 12*log2(MX(:,mask,13)/440), MX(:,mask,48), 'b', 'filled');
grid on; box on;
xlabel('semitones wrt 440');ylabel('max(a2/a1 smooth) during trans.');

hold on;
plot( 12*log2(F1MEAN(mask)/440),median(MX(:,mask, 48),1,'omitnan') ,  '-ro');



figure();
boxplot(MX(:,mask,48), 'Positions', 12*log2(F1MEAN(mask)/440), 'Widths',1);
box on; grid on;
xlabel({'xticks sample pipe','spacing semitones wrt 440'});
ylabel('max(a2/a1 smooth) during trans.'); 
ylim([0 50]);



%% Transient integrated spectral area between II and I
figure();
scatter( 12*log2(MX(:,:,13)/440),    MX(:,:,38).*MX(:,:,13)   ,'b' , 'filled');
hold on;
plot( 12*log2(F1MEAN/440), median( MX(:,:,38).*MX(:,:,13), 1, 'omitnan')  , '-ro');
ylabel('T1 * trapez(dt, a_2 - a_1), (t20_p, t80_p) '); xlabel('tessitura semitones');
xlabel('tessitura');
grid on; box on;
title('Integrated area comprised between a2 and a1 in the interval (t^{20}_p, t^{80}_p) normalized by oscill. period');
ylim([-5 40]);

%% Transient integrated (spectral) in area between III and I
figure();
scatter( 12*log2(MX(:,:,13)/440),    MX(:,:,39).*MX(:,:,13)   ,'b' , 'filled');
% scatter( 12*log2(MX(:,:,13)/440), MX(:,:,13)   ,'b' , 'filled');
hold on;
plot( 12*log2(F1MEAN/440), median( MX(:,:,39).*MX(:,:,13), 1, 'omitnan')  , '-ro');

ylabel('T1 * trapez(dt, a_3 - a_1), (t20_p, t80_p) '); xlabel('tessitura semitones');

grid on; box on;
title('Integrated area comprised between a3 and a1 in the interval (t^{20}_p, t^{80}_p) normalized by oscill. period');
ylim([-5 20]);







%%

fs = 10e3;
dt = 1/fs;
t = 0:dt: 0.75;
f = 1;
w = 2*pi*f;
x = 1*sin(w*t +pi/2) + 1;
xx = ones(size(x));

figure();
plot(t, x); grid on; hold on;
plot(t, xx);
tot = x-xx;

int = trapz(dt, tot)




%% Fig6, CFA Ernoult2016

Uj_at_a2max = sqrt(MX(:,:,46)*2/1.2);
rho = 1.2;
c = 340;
theta_a2max = Uj_at_a2max./(MX(:,:,9).*MX(:,:,13));

a2_nondim = MX(:,:,47)./(rho*c*Uj_at_a2max);

figure();
% scatter( 12*log2(MX(:,:,13)/440),a2_nondim);
scatter(theta_a2max, 1e3*a2_nondim,'filled');
xlabel('theta at t=a2max');
ylabel('a2max / rho c Uj at a2max');
grid on; box on;
xlim([0 11]);



%% a2max over a2 target
figure();
scatter( 12*log2(MX(:,:,13)/440), MX(:,:,45), 'b', 'filled');
xlabel('freq in semitones');
ylabel('a2 max over a2 target');
grid on; box on;ylim([0 5]);

%% Pmouth as per massage to equations 9-10-11 and alpha_vc = 1;


Pm =  MX(:,:,21) - (MX(:,:,5)./MX(:,:,6)).^2.*(MX(:,:,20) - MX(:,:,21));
Pm =  MX(:,:,21) + (MX(:,:,5)./MX(:,:,6)).^2.*(MX(:,:,21) - MX(:,:,20));

figure();

scatter(  12*log2(MX(:,:,13)/440), log10(Pm) , 'filled'); box on; grid on;

xlabel('$12log_2(f_1/440)$','interpreter','latex'); box on; grid on;
ylabel('P_{m} [Pa]','interpreter','latex');


%% t20rad - t20 foot [ACOUSTIC DELAY][OK]

figure();
scatter( 12*log2(MX(:,:,13)/440), 1e3*(MX(:,:,35)-MX(:,:,34) ), 'b', 'filled');
ylabel('$t^{20}_{rad}$ - t$^{20}_{foot}$ (lin)[ms]','interpreter','latex');
xlabel('$12log_2(f_1/440)$','interpreter','latex'); box on; grid on;
title('Acoustic delay','interpreter','latex');ylim([0 100])

%% (t20rad-t20foot)/T1 [ACOUSTIC DELAY][OK]

figure();
scatter( 12*log2(MX(:,:,13)/440), (MX(:,:,35)-MX(:,:,34) ).*MX(:,:,13), 'b', 'filled');
ylabel('$\frac{t^{20}_{rad} - t^{20}_{foot}}{T_1}$','interpreter','latex', 'Rotation',0);
xlabel('$12log_2(f_1/440)$','interpreter','latex'); box on; grid on;
title('Acoustic delay','interpreter','latex');%ylim([0 ])

%% t20_foot [HYDRODYNAMIC DELAY][OK]

figure();
scatter( 12*log2(MX(:,:,13)/440),1e3*(MX(:,:,34) ), 'b', 'filled');
ylabel('t$^{20}_{foot}$ (lin)[ms]','interpreter','latex');
xlabel('$12log_2(f_1/440)$','interpreter','latex'); box on; grid on;
title('Hydrodynamic delay','interpreter','latex');ylim([0 10])

%% t20_foot normalized by tau_f =  [HYDRODYNAMIC DELAY]

figure();
scatter( 12*log2(MX(:,:,13)/440),( c * MX(:,:,34).*MX(:,:,5)./MX(:,:,2) ), 'b', 'filled');
ylabel('$ \frac{t^{20}_f}{\tau_f} = \frac{c_o \times t^{20}_{f} \mathcal{S}_{in}}{V_f}$','interpreter','latex','Rotation',0);
xlabel('$12\times log_2(f_1/440)$','interpreter','latex'); box on; grid on;
title('Hydrodynamic delay','interpreter','latex');



%%
figure();
scatter(MX(:,:,29).*MX(:,:,13), (MX(:,:,43)), 'b','filled' );

%% Sj / Sin [OK]

figure();
scatter( 12*log2(MX(:,:,13)/440), log10(MX(:,:,6)./MX(:,:,5) ), 'b', 'filled');
ylabel('$\mathcal{S}_j / \mathcal{S}_{in}$ (log10)');
xlabel('$12log_2(f_1/440)$'); box on; grid on;


%% a2max/a1(@a2max) (mouth rad) [OK]
figure();
scatter(12*log2(MX(:,:,13)/440), (MX(:,:,43)), 'b','filled' );
xlabel({'Sample pipe xticks','12log2(f1/440) spacing'},'interpreter','latex');
ylabel('a2max/a1(@a2max) (log$_{10}$)','interpreter','latex'); 
box on; grid on;

% f1lump

figure();
funh = boxplot( log10(abs(MX(:,:,43))) , 'Positions',12*log2(F1MEAN/440), 'Widths',1);
% xlabel('Sample pipe','interpreter','latex');
% ylabel('a2max/a1(@a2max) (log$_{10}$)','interpreter','latex');



%% a2max/a1target (mouth rad) [OK]

figure();
scatter(12*log2(MX(:,:,13)/440), (MX(:,:,44)), 'b','filled' );
xlabel('$12*log_2(f_1/f_{ref})$','interpreter','latex');
ylabel('a2max/a1target (lin)','interpreter','latex'); box on; grid on;

figure();boxplot( (MX(:,:,44)), 'Positions',12*log2(F1MEAN/440), 'Widths',1);
xlabel({'sample pipe xticks','12 log2(f1/440) spacing'},'interpreter','latex');
ylabel('a2max/a1target (lin)','interpreter','latex');




%% PRTfoot/T1 (w.r.t tessitura) [T1 for normalization??? it's hydrodynamic]
figure();
scatter(12*log2(MX(:,:,13)/440), MX(:,:,29).*MX(:,:,13), 'b','filled' ); 
xlabel('$12\times log_2(f_1/f_{ref})$','interpreter','latex');
ylabel('$\frac{PRT_{foot}}{T_1}$','interpreter','latex','Rotation',0);
box on; grid on;

%% PRTfoot/T1 log10 (w.r.t tessitura) [T1 for normalization??? it's hydrodynamic]
figure();
scatter(12*log2(MX(:,:,13)/440), log10(MX(:,:,29).*MX(:,:,13)), 'b','filled' ); 
xlabel('$12\times log_2(f_1/f_{ref})$','interpreter','latex');
ylabel('$\frac{PRT_{foot}}{T_1}$ (log10)','interpreter','latex','Rotation',0);
box on; grid on;

%% PRTrad/T1 (w.r.t tessitura) [OK]
figure();
scatter(12*log2(MX(:,:,13)/440), MX(:,:,30).*MX(:,:,13), 'b','filled' ); 
xlabel('$12 \times log_2(f_1/f_{ref})$','interpreter','latex');
ylabel('$\frac{PRT_{rad}}{T_1}$','interpreter','latex','Rotation',0);
box on; grid on;

%% beta vs PRTfoot [OK]
figure();
scatter( 1./MX(:,:,26), MX(:,:,29)  , 'b','filled');
xlabel('$\beta^{-1}$', 'interpreter','latex');
ylabel('PRT$_{foot}$ [s]', 'Interpreter','latex');
box on; grid on;
xlim([0 2.2e-3]);
ylim([0 4.5e-3]);

%% beta vs exp(nu) [NO]
figure();
mask = [1:10,13:22];
mask = 1:22;
% scatter( MX(:,mask,26), (MX(:,mask,27)), 'b', 'filled');
% scatter( log10(MX(:,mask,26))/log10(20), log(MX(:,mask,27)), 'b', 'filled');
scatter( (MX(:,mask,26)), exp(MX(:,mask,27)), 'b', 'filled');
xlabel('$\beta$', 'interpreter','latex');
ylabel('$\nu$', 'Interpreter','latex');

%% nu vs beta^(0.56) [OK]
figure();
mask = [1:10,13:22];
mask = 1:22;
% scatter( MX(:,mask,26), (MX(:,mask,27)), 'b', 'filled');
% scatter( log10(MX(:,mask,26))/log10(20), log(MX(:,mask,27)), 'b', 'filled');
scatter( MX(:,mask,27), MX(:,mask,26).^(0.56), 'b', 'filled');
ylabel('$\beta^{0.56}$', 'interpreter','latex');
xlabel('$\nu$', 'Interpreter','latex'); grid on; box on;

%% Ptarg differences [OK]
figure();

scatter(12*log2(MX(:,:,13)/440), abs(MX(:,:,19)-MX(:,:,20)) ,'blue','filled');
figure();
scatter(12*log2(MX(:,:,13)/440), abs(MX(:,:,20)-MX(:,:,21)) ,'red','filled');
figure();
scatter(12*log2(MX(:,:,13)/440), abs(MX(:,:,21)-MX(:,:,22)), 'green','filled');

xlabel('$12 \times log_2(f_1/f_{ref})$','interpreter','latex');
ylabel('Target pressure differences','interpreter','latex');
% legend();



%%  Beta x T1 vs freq [OK]

mask = [1:22]; % low: sample pipes 2 and 3; high, sample pipes 11 and 12

figure();
% scatter( 12*log2(MX(:,mask,13)/440) , log(MX(:,mask,26)./MX(:,mask,13)) ,'b','filled');ylabel('beta x T_1 (log)');xlabel('12log_2(f_1/440)'); box on; grid on; 
scatter( 12*log2(MX(:,mask,13)/440) , (MX(:,mask,26)./MX(:,mask,13)) ,'b','filled');
ylabel('$\beta \times  T_1$','interpreter','latex');
xlabel('$12log_2(f_1/f_{ref})$','interpreter','latex'); 
box on; grid on; 

%% Beta x T1 vs freq BOX plot [OK]

figure();
boxplot(MX(:,:,26)./MX(:,:,13), 'Positions',12*log2(F1MEAN/440), 'Widths',1);
xlabel({'Sample pipe xticks','12log_2() fspacing'},'interpreter','latex');
ylabel('$\beta \times T_1$','interprete','latex');

%% nu [OK]
figure();
scatter( 12*log2(MX(:,:,13)/440) , log10(MX(:,:,27)) , 'b','filled'); ylabel('nu log10');

figure(); boxplot( log10(MX(:,:,27)) , 'Positions', 12*log2(F1MEAN/440), 'Widths', 1); ylabel('nu log10');

%% BETA vs Anything
namevarsall = {'Lp','Vf','PWdth','TNHD','INLET','SJET','h','H','Wm','Dp',...
    'Vgrv','Qfact','f1','theta','Qpall2grv','Qgr2ft','Qjet','Rmplssg',...
    'P{o} pall','P{o}grv','P{o}foot','P{o}rad',...
    'P{o} grv-pall','P{o}foot-grv','P{o}rad-foot',...
    'Beta','nu',...
    'PRTgrv','PRTfoot','PRTrad',...
    'PRTfoot/grv','PRTrad/foot',...
    't{20}grv','t{20}foot','t{20}rad',...
    't20(foot-grv)','t20(rad-foot)',...
    'Area1','Area2','Sin/Spal','Sjet/Sin','MaxKeyVel','a2max/a1(@a2max)','a2max/a1target'};
figure();
for idx = [1:44]
   scatter( log10(MX(:,:,26)), log10(MX(:,:,idx)) , 'b','filled');
   xlabel('beta (log10)'); ylabel([namevarsall{idx},' (log10)']);%ylabel(num2str(idx));
   box on; drawnow();
   pause();
end

%% BETA*T1 versus NU [NO]
figure();
scatter(log10(MX(:,:,26)./MX(:,:,13)), log10(MX(:,:,27)), 'b','filled');
xlabel('beta*T1'); ylabel('nu');
%% BETA versus Qj [OK]
figure();
scatter( log10(MX(:,:,26)), log10(MX(:,:,17)), 'b','filled');
xlabel('beta log10'); ylabel('Qj log10');
%% BETA*T1 versus Qj [NO]
figure();
scatter( log10(MX(:,:,26)./MX(:,:,13)), log10(MX(:,:,17)), 'b','filled');
xlabel('beta*T1'); ylabel('Qj');
%% Beta*Vf,Qj
figure();
scatter( log10( MX(:,:,26) ), log10(MX(:,:,17)) ,'b','filled');
%% NU vs Anything (pending to check) [NO]
figure();
for idx = [43:44]
   scatter( log10(MX(:,:,idx)) ,log10(MX(:,:,27)), 'b','filled');
   try
    ylabel('nu log'); xlabel([namevarsall{idx},' (ln)']);%ylabel(num2str(idx));
   end
   box on; drawnow();
   pause();
end

%% PRTgrv * f1 vs NU [NO]
figure();
scatter( log(MX(:,:,30).*MX(:,:,13)), log(MX(:,:,27)), 'b', 'filled');
%% PRTgrv * f1 vs NU
figure();
scatter( log(MX(:,:,30).*MX(:,:,13)), log(MX(:,:,27)), 'b', 'filled');



%% MaxKeyVelocity vs Anything [NO]
figure(24); clf;
for idx = [1:44]
   scatter( log(MX(:,:,42)), log(abs(MX(:,:,idx))) , 'b','filled');
   try
    xlabel('VeloMaxAbs log'); ylabel([namevarsall{idx},' (ln)']);%ylabel(num2str(idx));
   end 
   box on; drawnow();
   pause();
end



%% BETA vs MaxKeyVel [NO]
figure();
scatter( MX(:,:,26), MX(:,:,42), 'b','filled');
xlabel('beta');ylabel('maxKeyVelo [m/s]');

%% NU vs MaxKeyVel [NO]
figure();
scatter( log(1.*MX(:,:,27)), log(1.*MX(:,:,42)), 'b','filled');
xlabel('nu');ylabel('maxKeyVelo [m/s]'); box on;




%% PRT rad / PRT foot
figure();
scatter( 12*log2(MX(:,:,13)/440), log10(MX(:,:,32) ), 'b', 'filled');
ylabel('PRTrad / PRT_{foot} (log10)');
xlabel('$12log_2(f_1/440)$');
%% PRT foot / PRT groove
figure();
scatter( 12*log2(MX(:,:,13)/440), log10(MX(:,:,31) ), 'b', 'filled');
ylabel('PRT$_{foot}$ / PRT$_{grv}$ (log10)');
xlabel('$12log_2(f_1/440)$');

%% t20rad - t20 groove [OK]

figure();
scatter( 12*log2(MX(:,:,13)/440), log10(MX(:,:,35)-MX(:,:,33) ), 'b', 'filled');
ylabel('t$^{20}_{rad}$ / t$^{20}_{grv}$ (log10)');
xlabel('$12log_2(f_1/440)$');

box on; grid on;% ylim([0 2.4e-3]);

%% t80 foot

figure();
scatter( 12*log2(MX(:,:,13)/440), 1e3*(MX(:,:,29)+MX(:,:,34) ), 'b', 'filled');
ylabel('$t^{80}_{foot}$ [ms]');
xlabel('$12log_2(f_1/440)$');



%% t80 rad

figure();
scatter( 12*log2(MX(:,:,13)/440), 1e3*(MX(:,:,30)+MX(:,:,35) ), 'b', 'filled');
ylabel('$t^{80}_{rad}$ [ms]');
xlabel('$12log_2(f_1/440)$');


%% VC factor groove [OK][But what does it mean?]

factor = (MX(:,:,5)./(MX(:,:,3)*PalletDepth)).*sqrt( (MX(:,:,20) - MX(:,:,21) )./(MX(:,:,19)-MX(:,:,20)));

figure();
scatter( MX(:,:,5)*1e6, factor, 'b', 'filled');box on;
xlabel('$S^{geo}_{in} \ [mm^2]$','interpreter','latex','fontsize',FSZ);
ylabel('$\Gamma_{grv}$','interpreter','latex','fontsize',FSZ);





%% #################################################################### 
% Ratio of effective areas as a function of geometric areas ratio
% #################################################################### 


FSZ = 22;
figure();

% S foot/pallet ====

datax = (MX(:,:,5))./(MX(:,:,3)*PalletDepth); % Sin/Sgroove GEOMETRIC section
datay = sqrt( (  MX(:,:,19)-MX(:,:,20)  )./(MX(:,:,20)-MX(:,:,21)) ); % Effective section (SS)

axeshandle(1) = subplot(121);
scatter( datax , datay , 'b', 'filled');box on;
xlabel('$\frac{S^{\ geo}_{foot}}{S^{\ geo}_{pall}}$','interpreter','latex','fontsize',FSZ+8);
ylabel('$\frac{S^{\ eff}_{foot}}{S^{\ eff}_{pall}}$','interpreter','latex','fontsize',FSZ);
set(get(gca,'ylabel'),'Rotation',0);
grid on; ylim([0 2]);
set( get(gca,'XAxis'), 'FontSize', FSZ-8);
set( get(gca,'YAxis'), 'FontSize', FSZ-8);

% S jet/foot ====

dataone = MX(:,:,6)./MX(:,:,5); % GEOMETRIC section Sjet/Sin,foot
datatwo = sqrt( (MX(:,:,20)-MX(:,:,21)   )./( MX(:,:,21)  ));


axeshandle(2) = subplot(122);
scatter( dataone, datatwo, 'b', 'filled');box on;
hold on;
plot([0,1.4],[0,1.4],'--r');
xlabel('$\frac{S^{geo}_{jet}}{S^{geo}_{foot}}$','interpreter','latex','fontsize',FSZ);
ylabel('$\frac{S^{eff}_{jet}}{S^{eff}_{foot}}$','interpreter','latex','fontsize',FSZ);
set(get(gca,'ylabel'),'Rotation',0);
grid on; ylim([0 2]);
set( get(gca,'XAxis'), 'FontSize', FSZ-8);
set( get(gca,'YAxis'), 'FontSize', FSZ-8);

%%
% #################################################################### 
% GOOD: physical meaning and expected values
% #################################################################### 

%(VC factor foot)

FSZ = 22;
dataone = MX(:,:,6)./MX(:,:,5); % GEOMETRIC section Sjet/Sin,foot
datatwo = sqrt( (MX(:,:,20)-MX(:,:,21)   )./( MX(:,:,21)  ));

figure();
scatter( MX(:,:,5)*1e6, dataone./datatwo, 'b', 'filled');box on;
% scatter( dataone, dataone./datatwo, 'b', 'filled');box on;
% scatter( sqrt(MX(:,:,5)/pi)./MX(:,:,3), dataone./datatwo, 'b', 'filled');box on;
xlabel('$S^{geo}_{in} \ [mm^2]$','interpreter','latex','fontsize',FSZ);
ylabel('$\Gamma$ factor','interpreter','latex','fontsize',FSZ);
set(get(gca,'ylabel'),'Rotation',0);
% grid on; ylim([0 2]);
set( get(gca,'XAxis'), 'FontSize', FSZ);
set( get(gca,'YAxis'), 'FontSize', FSZ);
ylim([0,1]);
%%
figure(); 
% scatter(MX(:,:,3), MX(:,:,5), 'b', 'filled');
plot( (MX(:,:,3)*0.05./MX(:,:,5))' , 'o');

%%

figure();
scatter( MX(:,:,5), datatwo, 'b', 'filled');box on;
xlabel('$\frac{S^{geo}_{jet}}{S^{geo}_{foot}}$','interpreter','latex','fontsize',FSZ);
ylabel('$\frac{S^{eff}_{jet}}{S^{eff}_{foot}}$','interpreter','latex','fontsize',FSZ);
set(get(gca,'ylabel'),'Rotation',0);
grid on; ylim([0 2]);
set( get(gca,'XAxis'), 'FontSize', FSZ);

%%
XX = datax .* dataone;
YY = datay .* datatwo;

figure();

scatter( XX, YY, 'b', 'filled');box on;
% xlabel('$\frac{S^{geo}_{jet}}{S^{geo}_{foot}}$','interpreter','latex','fontsize',FSZ);
% ylabel('$\frac{S^{eff}_{jet}}{S^{eff}_{foot}}$','interpreter','latex','fontsize',FSZ);
set(get(gca,'ylabel'),'Rotation',0);
grid on; ylim([0 2]);
set( get(gca,'XAxis'), 'FontSize', FSZ-8);


%% 
% #################################################################### 
% AUTO CORRELATION OF VARIABLE VALUES 
% #################################################################### 
tmp = MXresh;

for idx = 1 : size(tmp,1)
    for jdx = 1 : size(tmp,2)
        if isnan(tmp(idx,jdx))
            tmp(idx,jdx) = 0;
        end
    end
end

DataToWorkWith = MXresh(:,maskpca);


tmp = (tmp-mean(tmp,1,'omitnan'))./var(tmp,'omitnan');

% tmp = log(abs(tmp)+1e-6);
tmp = exp(tmp);



tmpc = corrcoef(tmp);

tmpc(abs(tmpc)<0.0001) = 0;

figure();
imagesc(tmpc);
ax=gca;
ax.YDir = 'normal';
ax.XTick = 1:length(tmpc);
ax.XTickLabel = namevarsall;
ax.YTick = 1:length(tmpc);
ax.YTickLabel = namevarsall;
% colormap(bluewhitered); 

colorbar;


%%
% #################################################################### 
% PCA()
% #################################################################### 

% maskpca = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39];
% maskpca = [2,3,5,6,7,26,38];


DataToWorkWith = MXresh(:,maskpca);


wei = (DataToWorkWith-mean(DataToWorkWith,1,'omitnan'))./var(DataToWorkWith,'omitnan');

wei = ones(length(maskpca),1)';

[coeff,score,latent,tsquared,explained] = pca(DataToWorkWith, 'VariableWeights',wei);
figure(3);clf; plot(explained,'o-');
LBL = append(namevarsall(maskpca));


figure(4);clf;
biplot(coeff(:,1:3),'scores',score(:,1:3),'VarLabels',LBL);

%% 
% #################################################################### 
% BETA and NU with FREQUENCY
% #################################################################### 

figure(25); clf;
FSZ = 15;
scatter( freqlogax, log(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
hold on;
yyaxis right;
scatter( freqlogax, log(MX(:,:,27)), 'r', 'filled');box on;
% ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
legend('$\beta$','$\nu$','interpreter','latex','fontsize',FSZ);
yyaxis right;

ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
xlabel('$12log_2(f_1 / f_{ref})$','interpreter','latex','fontsize',FSZ);

%% 
% #################################################################### 
% BETA and NU with FREQUENCY
% #################################################################### 

figure(25); clf;
FSZ = 15;
scatter( freqlogax, log(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
hold on;
yyaxis right;
scatter( freqlogax, log(MX(:,:,27)), 'r', 'filled');box on;
% ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
legend('$\beta$','$\nu$','interpreter','latex','fontsize',FSZ);
yyaxis right;

ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
xlabel('$12log_2(f_1 / f_{ref})$','interpreter','latex','fontsize',FSZ);

%% 

% #################################################################### 
% BETA and NU with I21 and I31
% #################################################################### 
figure(); clf;
FSZ = 15;
scatter( log(MX(:,:,38)), log(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
xlabel('I21');

figure();
scatter( log(MX(:,:,38)), log(MX(:,:,27)), 'r', 'filled');box on;
ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
xlabel('I21');

%% 
% #################################################################### 
% BETA and NU with PTARG GROOVE
% #################################################################### 


figure(27); clf;
FSZ = 15;
scatter( log(MX(:,:,20)), log(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
hold on;
yyaxis right;
scatter( log(MX(:,:,20)), log(MX(:,:,27)), 'r', 'filled');box on;
% ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
% legend('$\beta$','$\nu$','interpreter','latex','fontsize',FSZ);
yyaxis right;
ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
xlabel('$log PTARG groove$','interpreter','latex','fontsize',FSZ);



%% 
% #################################################################### 
% BETA and NU with PTARG FOOT
% #################################################################### 


figure(28); clf;
FSZ = 15;
scatter( log(MX(:,:,21)), log(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
hold on;
yyaxis right;
scatter( log(MX(:,:,21)), log(MX(:,:,27)), 'r', 'filled');box on;
% ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
% legend('$\beta$','$\nu$','interpreter','latex','fontsize',FSZ);
yyaxis right;
ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
xlabel('$log PTARG foot$','interpreter','latex','fontsize',FSZ);


%% 

% #################################################################### 
% PRESSURE DROPS
% #################################################################### 

FSZ = 14;

fili = [1:9,12,15:22];
fili = [1:22];
DROP1 = MX(:,fili,20)-MX(:,fili,19) ;
DROP2 = MX(:,fili,21)-MX(:,fili,20);
DROP3 = -MX(:,fili,21);%-MX(:,fili,19);

if 0
figure(24);clf;
sbh(1) = subplot(1,2,1);
sasa=pcolor(DROP1); colormap(inferno);colorbar; sasa.FaceColor = 'interp'; sasa.LineStyle = 'none';
title(sprintf('Pallet2Groove. Mean drop: %1.2f [Pa], std %1.3f [Pa]',mean(DROP1,[1,2],'omitnan'), std(DROP1,0,[1 2],'omitnan')   ));

sbh(2) = subplot(1,2,2);
soso=pcolor(DROP2); colormap(inferno);colorbar; soso.FaceColor = 'interp'; soso.LineStyle = 'none';
title(sprintf('Groove2Foot. Mean drop: %1.2f [Pa], std %1.3f [Pa]',mean(DROP2,[1,2],'omitnan'), std(DROP2,0,[1 2],'omitnan')   ));

linkaxes(sbh,'xy');
ylim([1 36])

figure(13); clf; 
sisi = pcolor(DROP3);colormap(inferno);colorbar;sisi.FaceColor='interp';sisi.LineStyle='none';ylim([1 38]);
title(sprintf('Pallet2Foot. Mean drop: %1.2f [Pa], std %1.3f [Pa]',mean(DROP3,[1,2],'omitnan'), std(DROP3,0,[1 2],'omitnan')   ));
end

figure(28);clf;
scatter(fili, DROP1,'b');
hold on;
scatter(fili, DROP2, 'r');
xlabel('Num pipe','interpreter','latex','fontsize',FSZ);
ylabel('Pressure drop [Pa]','interpreter','latex','fontsize',FSZ);
scatter(fili, DROP3, 'g');







%% 
% #################################################################### 
% VENA CONTRACTA Inlet-Jet (steady state values)
% #################################################################### 

FSZ = 16;



MXh   = MX(:,:,7);
MXH   = MX(:,:,8);
MXpal = MX(:,:,19);
MXgr  = MX(:,:,20);
MXft  = MX(:,:,21);

rat2 = (MXh.*MXH)./(MX(:,:,5)).*sqrt( (MXft)./(MXgr-MXft)  );

figure();
RatioToPlot = 1./rat2;
linidx = 1:min(size(RatioToPlot));
scatter(linidx, RatioToPlot, 'b','filled');%ylim([0 1]);
title('Vena contracta ratios: $VC_{jet}/VC_{inlet}$','interpreter','latex','fontsize',FSZ);
p = polyfit(linidx, mean(RatioToPlot,1,'omitnan'),1);
hold on;
plot(polyval(p,linidx),'--r');
xlabel('Num pipe','interpreter','latex','fontsize',FSZ);


%% 

% #################################################################### 
% VENA CONTRACTA palletbox-groove-foot
% #################################################################### 


PALLAREA = repmat(PALLAREA(maskpipes)',max(size(MXgr)),1);

rat1 = MX(:,:,5)./PALLAREA .*sqrt( (MXgr-MXft) ./ (MXpal-MXgr) );

figure();
RatioToPlot = 1./rat1;
linidx = 1:min(size(RatioToPlot));
scatter(linidx, RatioToPlot, 'b','filled');%ylim([0 1]);
title('Vena contracta ratios: $VC_{inlet}/VC_{groove-slot}$','interpreter','latex','fontsize',FSZ);
p = polyfit(linidx, mean(RatioToPlot,1,'omitnan'),1);
hold on;
plot(polyval(p,linidx),'--r');
xlabel('Num pipe','interpreter','latex','fontsize',FSZ);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FSZ = 17;
figure();


geoh(1) = subplot(3,4,1);
scatter( freqlogax, log(MX(:,:,1)), 'b', 'filled');box on;
ylabel('$Lp$','interpreter','latex','fontsize',FSZ);

geoh(2) = subplot(3,4,2);
scatter( freqlogax, log(MX(:,:,2)), 'b', 'filled');box on;
ylabel('$Vf$','interpreter','latex','fontsize',FSZ);

geoh(3) = subplot(3,4,3);
scatter( freqlogax, log(MX(:,:,3)), 'b', 'filled');box on;
ylabel('$PWdth$','interpreter','latex','fontsize',FSZ);

geoh(4) = subplot(3,4,4);
scatter( freqlogax, log(MX(:,:,4)), 'b', 'filled');box on;
ylabel('$TnHD$','interpreter','latex','fontsize',FSZ);

geoh(5) = subplot(3,4,5);
scatter( freqlogax, log(MX(:,:,5)), 'b', 'filled');box on;
ylabel('$Sin$','interpreter','latex','fontsize',FSZ);

geoh(6) = subplot(3,4,6);
scatter( freqlogax, log(MX(:,:,6)), 'b', 'filled');box on;
ylabel('$Sjet$','interpreter','latex','fontsize',FSZ);

geoh(7) = subplot(3,4,7);
scatter( freqlogax, log(MX(:,:,7)), 'b', 'filled');box on;
ylabel('$h$','interpreter','latex','fontsize',FSZ);

geoh(8) = subplot(3,4,8);
scatter( freqlogax, log(MX(:,:,8)), 'b', 'filled');box on;
ylabel('$H$','interpreter','latex','fontsize',FSZ);

geoh(9) = subplot(3,4,9);
scatter( freqlogax, log(MX(:,:,9)), 'b', 'filled');box on;
ylabel('$Wm$','interpreter','latex','fontsize',FSZ);

geoh(10) = subplot(3,4,10);
scatter( freqlogax, log(MX(:,:,10)), 'b', 'filled');box on;
ylabel('$Dp$','interpreter','latex','fontsize',FSZ);

% geoh(11) = subplot(3,4,11);
% scatter( freqlogax,log( MX(:,:,11)), 'b', 'filled');box on;
% ylabel('$Vgrv$','interpreter','latex','fontsize',FSZ);

geoh(12) = subplot(3,4,12);
scatter( freqlogax, log(MX(:,:,12)), 'b', 'filled');box on;
ylabel('$Qfact1$','interpreter','latex','fontsize',FSZ);


xlabel('$12log_2(f_1 / f_{ref})$','interpreter','latex','fontsize',FSZ);
linkaxes(geoh,'x');    


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures 2-3-4 (S-S) (pall-groove-foot)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FSZ = 17;

figure();
f1h(1)=subplot(8,2,1);
scatter( freqlogax, MX(:,:,15), 'b', 'filled');box on;
ylabel('$Q_{pall2gr}$','interpreter','latex','fontsize',FSZ); title('Steady-State');

f1h(2)=subplot(8,2,2); 
scatter( freqlogax, log(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
f1h(3)=subplot(8,2,3);
scatter( freqlogax, MX(:,:,16), 'b', 'filled');box on;
ylabel('$Q_{gr2ft}$','interpreter','latex','fontsize',FSZ);
f1h(4)=subplot(8,2,4);
scatter( freqlogax, MX(:,:,27), 'b', 'filled');box on;
ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
f1h(5)=subplot(8,2,5);
scatter( freqlogax, MX(:,:,19), 'b', 'filled');box on;
ylabel('$P^{o}_{pall}$','interpreter','latex','fontsize',FSZ);
f1h(6)=subplot(8,2,8);
scatter( freqlogax, MX(:,:,28), 'b', 'filled');box on;
ylabel('$PRT_{grv}$','interpreter','latex','fontsize',FSZ);
f1h(7)=subplot(8,2,7);
scatter( freqlogax, MX(:,:,20), 'b', 'filled');box on;
ylabel('$P^{o}_{grv}$','interpreter','latex','fontsize',FSZ);
f1h(8)=subplot(8,2,10);
scatter( freqlogax, MX(:,:,29), 'b', 'filled');box on;
ylabel('$PRT_{foot}$','interpreter','latex','fontsize',FSZ);
f1h(9)=subplot(8,2,9);
scatter( freqlogax, MX(:,:,21), 'b', 'filled');box on;
ylabel('$P^{o}_{foot}$','interpreter','latex','fontsize',FSZ);
f1h(10)=subplot(8,2,12);
scatter( freqlogax, MX(:,:,31), 'b', 'filled');box on;
ylabel('$PRT_{(foot/grv)}$','interpreter','latex','fontsize',FSZ);
f1h(11)=subplot(8,2,11);
scatter( freqlogax, MX(:,:,23), 'b', 'filled');box on;
ylabel('$P^{o}_{(grv-pall)}$','interpreter','latex','fontsize',FSZ);
f1h(12)=subplot(8,2,14);
scatter( freqlogax, MX(:,:,33), 'b', 'filled');box on;ylim([0 0.010]);
ylabel('$t^{20}_{grv}$','interpreter','latex','fontsize',FSZ);
f1h(13)=subplot(8,2,13);
scatter( freqlogax, MX(:,:,24), 'b', 'filled');box on;
ylabel('$P^{o}_{(ft-grv)}$','interpreter','latex','fontsize',FSZ);
% f1h(14)=subplot(8,2,14);
% scatter( freqlogax, MX(:,:,34), 'b', 'filled');box on;ylim([0 0.010]);
% ylabel('$t^{20}_{foot}$','interpreter','latex','fontsize',FSZ);
f1h(15)=subplot(8,2,16);
scatter( freqlogax, MX(:,:,36), 'b', 'filled');box on;ylim([0 0.003]);
ylabel('$t^{20}_{(ft-grv)}$','interpreter','latex','fontsize',FSZ);





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures 3-4-5 (S-S) (groove-foot-outside)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FSZ = 14;

figure();
f2h(1)=subplot(8,3,19);
scatter( freqlogax, MX(:,:,13), 'b', 'filled');box on;
ylabel('$f_1$','interpreter','latex','fontsize',FSZ);

f2h(2)=subplot(8,3,2);
scatter( freqlogax, MX(:,:,26), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient')

f2h(3)=subplot(8,3,6);
scatter( freqlogax, MX(:,:,33), 'b', 'filled');box on;
ylabel('$t^{20}_{grv}$','interpreter','latex','fontsize',FSZ);ylim([0 10e-3]);

f2h(4)=subplot(8,3,1);
scatter( freqlogax, MX(:,:,14), 'b', 'filled');box on;ylim([0 10]);
ylabel('$\theta$','interpreter','latex','fontsize',FSZ);title('Steady-State');

f2h(5)=subplot(8,3,3);
scatter( freqlogax, MX(:,:,27), 'b', 'filled');box on;
ylabel('$\nu$','interpreter','latex','fontsize',FSZ);title('Transient')

f2h(6)=subplot(8,3,9);
scatter( freqlogax, MX(:,:,34), 'b', 'filled');box on;
ylabel('$t^{20}_{foot}$','interpreter','latex','fontsize',FSZ);ylim([0 10e-3]);

f2h(7)=subplot(8,3,22);
scatter( freqlogax, MX(:,:,17), 'b', 'filled');box on;
ylabel('$Q_{jet}$','interpreter','latex','fontsize',FSZ);

f2h(8)=subplot(8,3,5);
scatter( freqlogax, MX(:,:,28), 'b', 'filled');box on;
ylabel('$PRT_{grv}$','interpreter','latex','fontsize',FSZ);

f2h(9)=subplot(8,3,12);
scatter( freqlogax, MX(:,:,35), 'b', 'filled');box on;ylim([0 0.03])
ylabel('$t^{20}_{rad}$','interpreter','latex','fontsize',FSZ);

f2h(10)=subplot(8,3,4);
scatter( freqlogax, MX(:,:,20), 'b', 'filled');box on;
ylabel('$P^{o}_{grv}$','interpreter','latex','fontsize',FSZ);

f2h(11)=subplot(8,3,8);
scatter( freqlogax, MX(:,:,29), 'b', 'filled');box on;
ylabel('$PRT_{foot}$','interpreter','latex','fontsize',FSZ);

f2h(12)=subplot(8,3,15);
scatter( freqlogax, MX(:,:,36), 'b', 'filled');box on;ylim([0 0.003]);
ylabel('$t^{20}_{(ft-grv)}$','interpreter','latex','fontsize',FSZ);

f2h(13)=subplot(8,3,7);
scatter( freqlogax, MX(:,:,21), 'b', 'filled');box on;
ylabel('$P^{o}_{foot}$','interpreter','latex','fontsize',FSZ);

f2h(14)=subplot(8,3,11);
scatter( freqlogax, MX(:,:,30), 'b', 'filled');box on;
ylabel('$PRT_{rad}$','interpreter','latex','fontsize',FSZ);

f2h(15)=subplot(8,3,18);
scatter( freqlogax, MX(:,:,37), 'b', 'filled');box on;ylim([0 0.02]);
ylabel('$t^{20}_{(rad-ft)}$','interpreter','latex','fontsize',FSZ);

f2h(16)=subplot(8,3,10);
scatter( freqlogax, MX(:,:,22), 'b', 'filled');box on;
ylabel('$P^{o}_{rad}$','interpreter','latex','fontsize',FSZ);

f2h(17)=subplot(8,3,14);
scatter( freqlogax, MX(:,:,31), 'b', 'filled');box on;
ylabel('$PRT_{(ft/grv)}$','interpreter','latex','fontsize',FSZ);

f2h(18)=subplot(8,3,[20,23]);
scatter( freqlogax, log(MX(:,:,38)), 'b', 'filled');box on;ylim([-9 -1]);
ylabel('$I1 \ (log)$','interpreter','latex','fontsize',FSZ);

f2h(19)=subplot(8,3,13);
scatter( freqlogax, MX(:,:,24), 'b', 'filled');box on;
ylabel('$P^{o}_{(ft-grv)}$','interpreter','latex','fontsize',FSZ);

f2h(20)=subplot(8,3,17);
scatter( freqlogax, MX(:,:,32), 'b', 'filled');box on;
ylabel('$PRT_{(rad/ft)}$','interpreter','latex','fontsize',FSZ);

f2h(21)=subplot(8,3,[21,24]);
scatter( freqlogax, log(MX(:,:,39)), 'b', 'filled');box on;ylim([-9 -1]);
ylabel('$I2 \ (log)$','interpreter','latex','fontsize',FSZ);

f2h(22)=subplot(8,3,16);
scatter( freqlogax, MX(:,:,25), 'b', 'filled');box on;
ylabel('$P^{o}_{(rad-ft)}$','interpreter','latex','fontsize',FSZ);


linkaxes(f2h,'x');
xlim([-21 23]);

    
%% %%%% PLOT %%%%%
    
% BENOIT'S PLOT (not the correct layers any more)
FSZ = 16;

figure(5); clf;
Bh(1) = nexttile([1 12]);
scatter( freqlogax, MX(:,:,33), 'b', 'filled');box on;
ylabel('$P^o_{foot}$','interpreter','latex','fontsize',FSZ);

Bh(2) = nexttile([1 12]);
scatter(freqlogax,  MX(:,:,6)./ MX(:,:,5), 'b', 'filled');box on;
ylabel('$S_{jet}/S_{in}$','interpreter','latex','fontsize',FSZ);

Bh(3) = nexttile([1 12]);
scatter(freqlogax, MX(:,:,16), 'b', 'filled');box on;
ylabel('$Q_{in}$','interpreter','latex','fontsize',FSZ);

Bh(4) = nexttile([1 12]);
scatter(freqlogax, MX(:,:,14), 'b', 'filled');box on;
ylabel('$\theta$','interpreter','latex','fontsize',FSZ);

Bh(5) = nexttile([1 12]);
scatter(freqlogax, MX(:,:,24), 'b', 'filled');box on;
ylabel('$GRV_{PRT}$','interpreter','latex','fontsize',FSZ);

Bh(6) = nexttile([1 12]);
scatter(freqlogax, MX(:,:,25), 'b', 'filled');box on;
ylabel('$FOOT_{PRT}$','interpreter','latex','fontsize',FSZ);

Bh(7) = nexttile([1 12]);
scatter(freqlogax, MX(:,:,23), 'b', 'filled');box on;
ylabel('$\nu$','interpreter','latex','fontsize',FSZ)

Bh(8) = nexttile([1 12]);
scatter(freqlogax, MX(:,:,22), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);

xlabel('$12log_2(f_1 / f_{ref})$','interpreter','latex','fontsize',FSZ);
linkaxes(Bh,'x');    
    
    
    
    
    
    
    
    
    
    
    
    
%%
    


scatter(axh(1),freqlogax, MX(:,:,1),'filled','b');
scatter(axh(2),freqlogax, MX(:,:,2),'filled','b');
scatter(axh(3),freqlogax, MX(:,:,3),'filled','b');

scatter(axh(4),freqlogax, MX(:,:,4),'filled','b');

scatter(axh(5),freqlogax, MX(:,:,5),'filled','b');
scatter(axh(6),freqlogax, MX(:,:,6),'filled','b');

scatter(axh(7),freqlogax, MX(:,:,7),'filled','b');ylim(axh(7),[0.5 1]);
scatter(axh(8),freqlogax, MX(:,:,8),'filled','b');ylim(axh(8),[0 1]);
scatter(axh(9),freqlogax, MX(:,:,9),'filled','b'); ylim(axh(9), [0 0.25]);
scatter(axh(10),freqlogax, MX(:,:,10),'filled','b');%ylim(axh(10),[0 3]);

scatter(axh(11),freqlogax, MX(:,:,11),'filled','b');%ylim(axh(11),[-0 0.4]);
scatter(axh(12),freqlogax, MX(:,:,12),'filled','b');ylim(axh(12),[1.3 2.1]);
scatter(axh(13),freqlogax, MX(:,:,13),'filled','b');ylim(axh(13),[0 30]);
scatter(axh(14),freqlogax, MX(:,:,14),'filled','b');
scatter(axh(15),freqlogax, MX(:,:,15),'filled','b');

scatter(axh(16),freqlogax,  MX(:,:,16),'filled','b');
scatter(axh(17),freqlogax,  MX(:,:,17),'filled','b');
scatter(axh(18),freqlogax,  MX(:,:,18),'filled','b');

scatter(axh(19),freqlogax, MX(:,:,19),'filled','b');

linkaxes(axh(:),'x');
xlim([0.99*min(freqlogax),1.05*max(freqlogax)]);
grid(axh, 'on');box(axh, 'on');
drawnow();


xlim([-20 25]);


%% 

% #################################################################### 
% PRT alones 
% #################################################################### 

figure();
nexttile([1,4]);   
scatter(freqlogax, PRTMX(:,:,1),'filled','b');%ylim([0.5 1]);
ylabel('groove'); ylim([0 0.005]);
nexttile([1,4]);   
scatter(freqlogax, PRTMX(:,:,2),'filled','b');%ylim([0.5 1]);
ylabel('foot');ylim([0 0.005]);
nexttile([1,4]);   
scatter(freqlogax, PRTMX(:,:,3),'filled','b');%ylim([0.5 1]);
ylabel('pipe');ylim([0 0.200]);
hold on;
plot(freqlogax, PRTpipeMedian,'ro-');


%% 

% #################################################################### 
% Q-factor vs BETA
% #################################################################### 

figure(12);clf; hold on;
% Q-factor vs Beta
try
for idx = 1 : length(QFAC1)
    scatter(QFAC1(idx), log(MX(:, idx, 5)), 'b','filled');
end
end
xlabel('$Q-$factor [s$^{-1}$]','interpreter','latex','fontsize',14);
ylabel('$\beta$ fitted (log) [s$^{-1}$]','interpreter','latex','fontsize',14);
hold on

BM = MX(:, :, 5);
BM = mean(BM,1, 'omitnan');
p = polyfit(QFAC1(find(BM)), log(BM), 1);
betalin = polyval(p, QFAC1(find(BM)) );
plot(QFAC1(find(betalin)), (betalin), '--r');

box on;
%% 

% #################################################################### 
% NU versus Sj/Sin
% #################################################################### 

figure(13); clf; hold on;

palletww_mask   = PW(maskpipes);
palletarea_mask = PALLAREA(maskpipes);
Vf_mask         = VF(maskpipes);
hmask           = h(maskpipes);
Wmmask          = WM(maskpipes);
Inletmask       = INLET(maskpipes);
Sjetmask        = SJET(maskpipes);

NU = MX(:, :, 27);
ratios = Sjetmask./Inletmask;
try
for idx = 1 : size(NU,2)
scatter(  log(Vf_mask(idx)), log(NU(:,idx)), 'b','filled'); % YES good
% scatter( log( ratios(idx) ), log((NU(:,idx))), 'b','filled'); 
% scatter( log( Inletmask(idx)./palletarea_mask(idx) ), log((NU(:,idx))),'b','filled'); 
% scatter( log( Sjetmask(idx)./palletarea_mask(idx) ), log((NU(:,idx))),'b','filled'); 
% scatter( log( Inletmask(idx) ), log((NU(:,idx))),'b','filled'); 

end
end


hold on;

NM = mean(NU,1, 'omitnan');
p = polyfit(log(ratios), log(NM), 1);
nulin = polyval(p, log(ratios) );
% plot(log(ratios), (nulin), '--r')


FSZ = 14;
xlabel('$S_{jet}/S_{inlet}$ (log)','interpreter','latex','fontsize',FSZ);
ylabel('$\nu$ (ln)','interpreter','latex','fontsize',FSZ);
box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% maskpca = [16,17,23,24,25,28,29,30,33,34,35,38,39];
maskpca = [16,17,19:22,28,29,30,33,34,35,38,39];
% maskpca = [16,17,24,29,30,33,34,35,38,39];

namevarsall = {'Lp','Vf','PWdth','TNHD','INLET','SJET','h','H','Wm','Dp',...
    'Vgrv','Qfact','f1','theta','Qpall2grv','Qgr2ft','Qjet','Rmplssg',...
    'P{o} pall','P{o}grv','P{o}foot','P{o}rad',...
    'P{o} grv-pall','P{o}foot-grv','P{o}rad-foot',...
    'Beta','nu',...
    'PRTgrv','PRTfoot','PRTrad',...
    'PRTfoot/grv','PRTrad/foot',...
    't{20}grv','t{20}foot','t{20}rad',...
    't20(foot-grv)','t20(rad-foot)',...
    'Area1','Area2','Sin/Spal','Sjet/Sin','MaxKeyVel'};

namevars2 = ['Lp','Vf','PWdth','TNHD','INLET','SJET','h','H','Wm','Dp',...
    'Vgrv','Qfact','f1','theta','Qpall2grv','Qgr2ft','Qjet','Rmplssg',...
    'PpallTarg','PgrvTarg','PfootTarg','PradTarg',...
    'Pgrv-Ppall(TARG)','Pfoot-Pgrv(TARG)','Prad-Pfoot(TARG)',...
    'PRTgrv','PRTfoot','PRTrad',...
    'PRTfoot/grv','PRTrad/foot',...
    't20grv','t20foot','t20rad',...
    't20(foot-grv)','t20(rad-foot)',...
    'Area1','Area2','Sin/Spal','Sjet/Sin'];

namevars3 = ["Lp","Vf","PWdth","TNHD","INLET","SJET","h","H","Wm","Dp",...
    "Vgrv","Qfact","f1","theta","Qpall2grv","Qgr2ft","Qjet","Rmplssg",...
    "PpallTarg","PgrvTarg","PfootTarg","PradTarg",...
    "Pgrv-Ppall(TARG)","Pfoot-Pgrv(TARG)","Prad-Pfoot(TARG)",...
    "PRTgrv","PRTfoot","PRTrad",...
    "PRTfoot/grv","PRTrad/foot",...
    "t20grv","t20foot","t20rad",...
    "t20(foot-grv)","t20(rad-foot)",...
    "Area1","Area2","Sin/Spal","Sjet/Sin"];

namevars4 = {"Lp","Vf","PWdth","TNHD","INLET","SJET","h","H","Wm","Dp",...
    "Vgrv","Qfact","f1","theta","Qpall2grv","Qgr2ft","Qjet","Rmplssg",...
    "P^{o}pall","P^{o}grv","P^{o}foot","P^{o}rad",...
    "P^{o}grv-Ppall","P^{o}foot-Pgrv","P^{o}rad-foot",...
    "PRTgrv","PRTfoot","PRTrad",...
    "PRTfoot/grv","PRTrad/foot",...
    "t20grv","t20foot","t20rad",...
    "t20(foot-grv)","t20(rad-foot)",...
    "Area1","Area2","Sin/Spal","Sjet/Sin"};



% With respect to beta
idxanalys = 26;
selmask   = [1:14]; % Only geometrical parameters 
namesmask = namevars4(selmask);
if 1

%     [idx, scores] = fsrftest(MXresh(:,[1:10,12:25,28:39]), MXresh(:,idxanalys));
    
    [idx, scores] = fsrftest(MXresh(:,selmask), MXresh(:,idxanalys));
    figure(1);clf;
    bar(scores(idx));
    set(gca, 'xtick',[1:length(idx)]);
    set(gca, 'xticklabels', namesmask(idx));
    title(sprintf('For: %s',namevars4{idxanalys}));
% else
    mdl = fsrnca(MXresh(:,selmask), MXresh(:,idxanalys));
    figure(2);clf;
    plot(mdl.FeatureWeights,'ro');
    grid on;
    xlabel('Feature index');
    ylabel('Feature weight');
    IDX = find( mdl.FeatureWeights > (mean(mdl.FeatureWeights)) );
    hold on;
    for jdx =1:length(IDX)
        text(IDX(jdx),mdl.FeatureWeights(IDX(jdx)),namesmask{IDX(jdx)});
    end
    xlim([0,length(mdl.FeatureWeights)+1]);
    
end 
% % % [idx, scores] = fsrmrmr(MXresh(:,[1:10,12:25,28:39]), MXresh(:,26));