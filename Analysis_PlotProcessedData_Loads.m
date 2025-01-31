% ============= LOADS =========================

thepath  = './processed/';
filecard = 'A*PROCESSED.mat';
files    = dir([thepath,filecard]);


% Resonator length: (56)
load('./processed/Geometry/Lp_m.mat'); 
LP = Lp_m; clear Lp_m

% Foot volume: (56)
load('./processed/Geometry/Vf_m3.mat'); 
VF = Vf_m3; clear Vf_m3

% Width of pallet window SLOT (56)
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

% Actual pallet valve dimensions (56 values)
load('./PalletValveDimensions.mat');
PVG = PalletValveGeometry; clear PalletValveGeometry;

% ============= PARAMETERS =========================

rho   = 1.2;
c     = 340;
maskpipes  = [3,4,5,6,7,9,10,11,13,15,17,19,24,25,27, 29,32, 34,37,39,41,44]';



PalletWinDepth = 1e-3*129.8; %  [m]
PALLwinAREA    = PalletWinDepth * PW;  % <<Slot>> area covered by the valve [m^2]
VGROOVE     = PW *0.51*0.05;

PalletValveWidth = PVG(:,2);
PalletValveVerticalMaxDisplacement = PVG(:,4);
PalletLength = 0.157;
PalletValveStrokeArea = (PalletValveWidth + PalletLength).*PalletValveVerticalMaxDisplacement ;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefs and params 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vena contracta(s):
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

n = 1; % mode / Method Blanc2010
chi       = 3e-5 * sqrt(n*F1MEAN)./RADSmask;
OM        = n*pi*c./(LPmask + 1.2*RADSmask);
part1     = OM/(2*c);
part2     = (LPmask + 1.2*RADSmask)./(chi.*LPmask + (OM.^2.*RADSmask.^2)/(2*c^2) )  ;
QFAC1     = part1.*part2;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MX            = nan*ones(NumTransMax, length(files), 58);
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
    onePalletValveStrokeArea = PalletValveStrokeArea(numpipe);
    
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

    MX(:, idx, 50)   = onePalletValveStrokeArea*ones(NumTransMax,1);
    
    
    
    %%%%%%%%%% S-S ASPECTS %%%%%%%%%% 
    
    % S1/S2 factor, Spall/SJET(Pf,Pgrv) Remplissage factor:
    Remplissage = real(sqrt(Pfoot_targ ./ (Pgroove_targ-Pfoot_targ ) ) ); 
    Qpall2grv   = VCpallet * PALLwinAREA(numpipe) * real(sqrt( (PpalletB_targ - Pgroove_targ) *2/rho));
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
    PRT20groove = PRT20groove(:);PRT20foot=PRT20foot(:);PRT20pipe=PRT20pipe(:);
    
    MX(find(betafit)  , idx, 26)  = betafit(:); % Foot pressure transient slope
    MX(find(nufit)    , idx, 27)  = nufit(:);   % Foot pressure transient onset sharpness
    
    
    MX(find(PRT20groove), idx,28) = PRT20groove;
    MX(find(PRT20foot)  ,idx,29)  = PRT20foot;
    MX(find(PRT20pipe)  ,idx,30)  = PRT20pipe;
    
    MX(find(PRT20foot)     ,idx, 31) = PRT20foot./PRT20groove;
    MX(find(PRT20pipe)     ,idx, 32) = PRT20pipe./PRT20foot;
    
    MX(find(t20groove), idx, 33) = t20groove(t20groove~=0);
    MX(find(t20foot),   idx, 34) = t20foot(t20foot~=0);
    MX(find(t20mouth),  idx, 35) = t20mouth(t20mouth~=0);
    
    MX(find(t20foot-t20groove), idx, 36) = (t20foot-t20groove);
    MX(find(t20mouth-t20foot),  idx, 37) = (t20mouth-t20foot);
    
    MX(find(betafit)            , idx, 38)   = Area1(:);
    MX(find(betafit)            , idx, 39)   = Area2(:);
    
    
    MX(:           , idx, 40)   = ones(NumTransMax,1)*(INLET(numpipe)/PALLwinAREA(numpipe));
    MX(:           , idx, 41)   = ones(NumTransMax,1)*(SJET(numpipe)/INLET(numpipe));
    
    MX( : , idx, 42)            = KeyVel(:,idx);

    MX(find(A2max_over_A1simult) ,idx,43)      = A2max_over_A1simult(find(A2max_over_A1simult));
    MX(find(A2max_over_A1target),idx,44)       = A2max_over_A1target(find(A2max_over_A1target));
    MX(find(A2max_over_A2target),idx,45)       = A2max_over_A2target(find(A2max_over_A2target));
    MX(find(pf_at_a2max),idx,46)               = pf_at_a2max(find(pf_at_a2max)) - 0;%+0*pm_at_a2max; % DeltaP (foot-mouth) /!\ Decide whether you keep the Pm'
    MX(find(a2max_vec),idx,47)                 = a2max_vec(find(a2max_vec));
    MX(find(max_a2_over_a1), idx, 48)          = max_a2_over_a1(find(max_a2_over_a1));
    MX(find(gofr2), idx, 49)                   = gofr2;


    MX(find(t10groove),  idx,51)  = t10groove;
    MX(find(PRT10groove),idx,52)  = PRT10groove;
    

    MX(find(t10foot),    idx,55)  = t10foot;
    MX(find(PRT10foot),  idx,56)  = PRT10foot;
    MX(find(t10mouth), idx, 57)   = t10mouth;
    MX(find(PRT10pipe), idx, 58)  = PRT10pipe;
    

   
    %%%%%%%%%% MARGINAL ASPECTS %%%%%%%%%% 
    PRTpipeMedian(idx) = median(PRT20pipe(PRT20pipe~=0),'omitnan');
    vecmeanfreqs(idx)  = mean(f1);
    freqref            = 440;
    freqlogax          = 12*log2(vecmeanfreqs/freqref);
    
    
end   
