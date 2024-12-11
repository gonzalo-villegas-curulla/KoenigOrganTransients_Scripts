try
    cd /run/media/gvc/ExtremeSSD/OrganPipe2023-2024/DataTransients/
end 

clc, clear;
addpath('./processed/');


 run Analysis_PlotProcessedData_Loads.m


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
% [50] Lateral stroke area of pallet valve at max opening (smaller than Slot Area = PWidth*Length_win_slot =MX(:,:,3)*0.129;


rho = 1.2;
co  = 340;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Geometry
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Vf vs Vgroove [OK][there's something interesting/intriguing...]

figure();
subplot(121);
scatter( 12*log2(MX(:,:,13)/440), log10(MX(:,:,11)./MX(:,:,2)), 'b', 'filled');
box on; grid on;
ylabel('$Vfoot / Vgroove \ (log_{10})$ [n.u.]');
xlabel('Tessitura [semitones]');

subplot(122);
scatter(  1e3*(MX(:,:,11)), 1e3*(MX(:,:,2)) , 'b', 'filled' );
xlabel('$V_{gr} \ [1e3 \times m^3]$');
ylabel('$V_f \ [1e3 \times m^3]$');
grid on; box on; xlim([0 0.45]);

%% Vf vs Sin vs Sj
figure();
scatter( (MX(:,:,2)./MX(:,:,3))  ,  (MX(:,:,6)), 'b', 'filled');
grid on; box on;
xlabel('$V_f$ / PalletValve-Slot-Width  $[m^2]$');
ylabel('$S_j \ [m^2]$');
ylim([0 7e-5]);

yyaxis right;
scatter((MX(:,:,2)./MX(:,:,3))  ,  (MX(:,:,5)), 'filled');
ylabel('$S_{in} \ [m^2]$');
ylim([0 7e-5]);


%% Sj / Sin [OK]

figure();
scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan'), median((MX(:,:,6)./MX(:,:,5)),1,'omitnan'), 'b', 'filled');
ylabel('$\mathcal{S}_j / \mathcal{S}_{in} \ (log_{10})$ [n.u.]', 'interpreter','latex');
xlabel('$12log_2(f_1/440)$','interpreter','latex'); box on; grid on; ylim([0 1.5]);




%% Sj / Sin */ Vf^0.33 [interactions]

figure();
scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan'), median(   (MX(:,:,6)./MX(:,:,5))./MX(:,:,2).^(1/3)   ,1,'omitnan'), 'b', 'filled');
ylabel('$\mathcal{S}_j / \mathcal{S}_{in} \ (log_{10})$ [n.u.]', 'interpreter','latex');
xlabel('$12log_2(f_1/440)$','interpreter','latex'); box on; grid on;

%% S_palletvalve, S_slot, S_tonehole, S_in, S_j 


figure(); hold on; grid on; box on;

scatter(freqlogax, 1e6*MX(:,:,50), 'b', 'filled');       % mean = 8.4 cm^2 % pallet valve stroke
scatter(freqlogax, 1e6*MX(:,:,3)*0.1298, 'r', 'filled'); % mean = 16.3 cm^2 % pallet slot window
scatter(freqlogax, 1e6*pi*(0.5*MX(:,:,4)).^2,'k','filled' ); % tone hole area
scatter(freqlogax, 1e6*MX(:,:,5), 'g', 'filled'); % foot inlet
scatter(freqlogax, 1e6*MX(:,:,6), 'm', 'filled'); % foot outlet
ylabel('mm^2'); xlabel('tessitura [semitones]');

%% Spallvalve / Sin by freqs

figure();
scatter(freqlogax, MX(:,:,50)./MX(:,:,5), 'b', 'filled');


%% | Pgrv - Ppall|/Ppall versus sections [NO]
figure();
% scatter( MX(:,:,50) , MX(:,:,20)); ylim([0 820]);
% scatter( MX(:,:,50) , abs(MX(:,:,20)-MX(:,:,19))./MX(:,:,19)    ); ylim([0 1]);%xlim([0 12e-4]);
% scatter( MX(:,:,3) , abs(MX(:,:,20)-MX(:,:,19))./MX(:,:,19)    ); ylim([0 1]);%xlim([0 12e-4]);
% scatter( (pi*(0.5*MX(:,:,4)).^2)./MX(:,:,50) , abs(MX(:,:,20)-MX(:,:,19))./MX(:,:,19)    ); ylim([0 1]);%xlim([0 12e-4]);
scatter( median((pi*(0.5*MX(:,:,4)).^2)./MX(:,:,50),1,'omitnan') , median(abs(MX(:,:,20)-MX(:,:,19))./MX(:,:,19),1,'omitnan') ,'b','filled'  );% ylim([0 1]);box on; grid on;
xlim([0 0.22]);
ylim([0 1]);

%% Pressure drops versus cross sections [NO]

figure(); % Pall2Groove

scatter( median((pi*(0.5*MX(:,:,4)).^2)./MX(:,:,50),1,'omitnan') , median(abs(MX(:,:,20)-MX(:,:,19))./MX(:,:,19),1,'omitnan') ,'b','filled'  );% ylim([0 1]);box on; grid on;


figure(); % Pall2Foot

% scatter( median(MX(:,:,6),1,'omitnan')  , median(  abs(MX(:,:,21)-MX(:,:,19))./MX(:,:,19),1,'omitnan')  , 'b', 'filled'); 
scatter( median(MX(:,:,6)./MX(:,:,50),1,'omitnan')  , median(  abs(MX(:,:,21)-MX(:,:,19))./MX(:,:,19),1,'omitnan')  , 'b', 'filled'); 

%% Pressure drops versus cross sections
figure();hold on;
    errorbar(   median((pi*(0.5*MX(:,:,4)).^2)./MX(:,:,50),1,'omitnan')  ,    median(abs(MX(:,:,19)-MX(:,:,20))./MX(:,:,19) , 1, 'omitnan'),...
        std( abs(MX(:,:,19)-MX(:,:,20))./MX(:,:,19), [], 1,'omitnan'),...
        'b','linestyle','none','marker','o');
    
    errorbar(   median(MX(:,:,5)./(pi*(0.5*MX(:,:,4)).^2),1,'omitnan')  , median(abs(MX(:,:,20)-MX(:,:,21))./MX(:,:,19) , 1, 'omitnan'),...
        std( abs(MX(:,:,20)-MX(:,:,21))./MX(:,:,19), [], 1,'omitnan'),...
        'r','linestyle','none','marker','^');
    
    errorbar( median(MX(:,:,6)./MX(:,:,5),1,'omitnan')  , median(abs(MX(:,:,21)-0)./MX(:,:,19) , 1, 'omitnan'),...
        std( abs(MX(:,:,21)-0         )./MX(:,:,19), [], 1,'omitnan'),...
        'g','linestyle','none','marker','x');

% xlim([-20 25]);ylim([0 1]); grid on; box on;
% ylabel('Pressure drops normalized by $P_{rsv}$','interpreter','latex');
xlabel('$...$','interpreter','latex');
legend('$\Delta P_{pall2gr}/P_{pall} = |P^{\oplus}_{gr}-P^{\oplus}_{pall}|/P^{\oplus}_{pall}$',...
    '$\Delta P_{gr2ft}\ /P_{pall} = |P^{\oplus}_{ft}-P^{\oplus}_{grv}|/P^{\oplus}_{pall}$',...
    '$\Delta P_{ft2m}\ /P_{pall}  = |\langle P^{\oplus}_{m}\rangle-P^{\oplus}_{ft}|/P^{\oplus}_{pall}$','interpreter','latex');


%%
figure();

scatter( median( (MX(:,:,5))./MX(:,:,50),1,'omitnan')  , median(  abs(MX(:,:,20))./MX(:,:,19),1,'omitnan')  , 'b', 'filled'); 


        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Steady-State analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Foot flow conservation ["appealing" result, but something is wrong]
      %   S_in                     SS target    SS target
      %   geom                      Preserv      Pfoot 
Qin = 1.0*MX(:,:,5) .* sqrt( 2/rho*( MX(:,:,19) - MX(:,:,21) ) );

       % Geom                      SS target   SS target
       % S_j                        P_foot      P_mouth
Qj  = 1.0*MX(:,:,6).*sqrt( 2/rho * ( MX(:,:,21) -      0     ) );

if 0
    figure();
    scatter( 12*log2(  MX(:,:,13   )/440) , abs(Qin-Qj)./Qin ,'b','filled');
    % scatter( 12*log2(  MX(:,:,13   )/440) , Qj./Qin ,'b','filled');
    grid on; xlabel('tessitura [semitones]');ylabel('Qj/Qin [m^3/s]'); box on;
    ylim([0 1]);
end

if 0
    figure();
    errorbar( 12*log2(F1MEAN/440) , median( abs(Qin-Qj)./Qin,1,'omitnan'), std( abs(Qin-Qj)./Qin,1, 'omitnan'),'k', 'marker','x','linestyle','none');
    grid on; ylim([0 1]);
    xlabel('$12\times log_2 (f_1 / 440 \ Hz)$','interpreter','latex');
    ylabel('$|Q_{in}-Q_j|/Q_{in}$','interpreter','latex');
end

% Qj/Qin wrt to Sj/Sin
figure();
errorbar(  median(MX(:,:,6)./MX(:,:,5),1,'omitnan') , median( abs(Qj)./Qin,1,'omitnan'), std( abs(Qj)./Qin,1, 'omitnan'),'k', 'marker','x','linestyle','none');
% errorbar(  median(12*log2(MX(:,:,13)/440),1,'omitnan') , median( abs(Qj)./Qin,1,'omitnan'), std( abs(Qj)./Qin,1, 'omitnan'),'k', 'marker','x','linestyle','none');
% plot(  median(MX(:,:,6)./MX(:,:,5),1,'omitnan') , median( abs(Qj)./Qin,1,'omitnan'), '-o');
grid on; ylim([0 1]);
xlabel('$ \mathcal{S}_j / \mathcal{S}_{in}$','interpreter','latex');
ylabel('$Q_j/Q_{in}$','interpreter','latex');
% xlim([0 1.4]);




%% Theta
figure();

scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan'), median(MX(:,:,14),1,'omitnan'), 'b', 'filled');
grid on; box on;
xlabel('Tessitura [semitones]');
ylabel('\theta = u_j/W_m f_1'); ylim([0 12]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target pressures study
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Ptarg differences [OK] all wrt to reservoir, in SS grv and rsv have the same pressure [OK, final]
% then we loose half to enter the foot and the other half is lost in the
% mouth outing: so the model is mouthpressure=0

% we can neglect the press drop between resrv and groove
% there are 2 significant press drops, of equivalent amgnitude, to enter
% foot and to exit the foot

% figure();
% scatter(12*log2(MX(:,:,13)/440), abs(MX(:,:,19)-MX(:,:,20))./MX(:,:,19) ,'blue','filled');hold on;
% scatter(12*log2(MX(:,:,13)/440), abs(MX(:,:,19)-MX(:,:,21))./MX(:,:,19) ,'red','filled');hold on;
% scatter(12*log2(MX(:,:,13)/440), abs(MX(:,:,21)-0 )./MX(:,:,19), 'green','filled');
% xlabel('$12 \times log_2(f_1/f_{ref})$','interpreter','latex');
% ylabel('Target pressure differences [n.u.]','interpreter','latex');
% title('Ptarg differences: blue=reserv2grv, red=grv2foot,green=foottomouth. ALL scaled by resrv  press.');
% grid on; box on;

%% Ptarg differences [OK] all wrt to reservoir, in SS grv and rsv have the same pressure [final]
figure();hold on;
if 0
    errorbar( 12*log2(F1MEAN/440), median(abs(MX(:,:,19)-MX(:,:,20))./MX(:,:,19) , 1, 'omitnan'),...
        std( abs(MX(:,:,19)-MX(:,:,20))./MX(:,:,19), [], 1,'omitnan'),...
        'k','linestyle','none','marker','o');
    errorbar( 12*log2(F1MEAN/440), median(abs(MX(:,:,20)-MX(:,:,21))./MX(:,:,19) , 1, 'omitnan'),...
        std( abs(MX(:,:,20)-MX(:,:,21))./MX(:,:,19), [], 1,'omitnan'),...
        'k','linestyle','none','marker','^');
    errorbar( 12*log2(F1MEAN/440), median(abs(MX(:,:,21)-0)./MX(:,:,19) , 1, 'omitnan'),...
        std( abs(MX(:,:,21)-0         )./MX(:,:,19), [], 1,'omitnan'),...
        'k','linestyle','none','marker','x');
else
    plot( 12*log2(F1MEAN/440), median(abs(MX(:,:,19)-MX(:,:,20))./MX(:,:,19) , 1, 'omitnan'),...       
        'k','linestyle','none','marker','o');
    plot( 12*log2(F1MEAN/440), median(abs(MX(:,:,20)-MX(:,:,21))./MX(:,:,19) , 1, 'omitnan'),...
        'k','linestyle','none','marker','^');
    plot( 12*log2(F1MEAN/440), median(abs(MX(:,:,21)-0)./MX(:,:,19) , 1, 'omitnan'),...
        'k','linestyle','none','marker','x');
end

xlim([-20 25]);ylim([0 1]); grid on; box on;
% ylabel('Pressure drops normalized by $P_{rsv}$','interpreter','latex');
xlabel('$12\times log_2 (f_1/440 \ Hz)$','interpreter','latex');
legend('$\widetilde{\Delta P}_{pall2grv} = |P^{\oplus}_{gr}-P^{\oplus}_{pall}|/P^{\oplus}_{pall}$',...
    '$\widetilde{\Delta P}_{grv2ft}\  = |P^{\oplus}_{ft}-P^{\oplus}_{grv}|/P^{\oplus}_{pall}$',...
    '$\widetilde{\Delta P}_{ft2m}\ \   = |\langle P^{\oplus}_{m}\rangle-P^{\oplus}_{ft}|/P^{\oplus}_{pall}$','interpreter','latex');

%% Pfoot target versus Sj/Sin [ok]

 % [Can we predict PfootTarg just by Sj and Sin?]
 
 Prsv  = 820.83;
 Prsv  = 1;
 
 % [1] DATA measured: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Sgeom_ratio_measured = median(MX(:,:,6)./MX(:,:,5),1,'omitnan');
 Qj_over_Qin    = median(  (MX(:,:,6).*sqrt(2/rho*(MX(:,:,21)-0)))  ./  (MX(:,:,5).*sqrt(2/rho*(MX(:,:,20)-MX(:,:,21)))) ,...
           1,'omitnan');
 
 
 FitFlowConserv = polyfit( Sgeom_ratio_measured , Qj_over_Qin ,1);
  
 A = FitFlowConserv(1); B = FitFlowConserv(2);
 Gamma    = (A * S_axes + B); 
 GammaInv = 1./Gamma;
 
 Sgeom_measured = median(MX(:,:,6)./MX(:,:,5),1,'omitnan');
 Seff_measured  = 0.5477*Sgeom_measured + 0.1398;
 
%  gamma_measured = median(MX(:,:,6)./MX(:,:,5).*sqrt(MX(:,:,21)./abs(MX(:,:,21)-MX(:,:,19))),1,'omitnan');
 gamma_measured = Qj_over_Qin;
 
 
 % [2] MODEL:  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 S_axes      = linspace(0,1.5, 1e3);
 Pfoot_model = Prsv./(1 + (GammaInv.*S_axes).^2);
 x = linspace(0,2.0,1e2);
 y = 1./(1+x.^2);
 
 y = 1./(1 + (x./(A*x+B)).^2);

 figure();
 % NO:
%  plot( 1./gamma_measured .* Sgeom_ratio_measured, ...
%      median(MX(:,:,21)./MX(:,:,19),1,'omitnan') , ...     
%      'xk');

% Poor
 plot( 1./(A * Sgeom_ratio_measured + B).*Sgeom_ratio_measured , ...
     median(MX(:,:,21)./MX(:,:,19),1,'omitnan') , ...     
     'xk');
 % Best
  plot( Sgeom_ratio_measured , ...
     median(MX(:,:,21)./MX(:,:,19),1,'omitnan') , ...     
     'xk');
 
 hold on;
 plot(x, y, '--k');
 ylim([0 1]);
 
 grid on; box on; xlim([0 1.5]);
 ylabel('$P_f^{\oplus}/P_{pall}^{\oplus}$ [n.u.]', 'interpreter','latex');
 xlabel('$\frac{\mathcal{S}_j}{ \mathcal{S}_{in}}$ [n.u.]', 'interpreter','latex');


%% uj/uin (no vena contracta) [OK]

figure();
scatter( 12*log2(MX(:,:,13)/440) , MX(:,:,18), 'b', 'filled');
grid on; box on;
ylim([0 1]);
xlabel('Tessitura [semitones]'); title('Remplissage');
 ylabel('$u_j / u_{in} = \sqrt{P_f^o} / \sqrt{ P_{gr}^o - P_f^o} $','interpreter','latex');

 % Additional hypothesis on why we don't see the flow consrevation:
 % the foot is a big volume with a small channel at the entrance and output
 % so we have a jet at the entrance dissipated by turbulence
 % and a second at the output

 %% Pgroove target versus (1) Stonehole/Spalletchoke or (2) Sin/Spalletchoke or (3) Sin/Stonehole
 

 figure();
 scatter(  median( MX(:,:,50)./MX(:,:,6),1,'omitnan'), median(MX(:,:,19),1,'omitnan') , 'b', 'filled');
 grid on; box on; ylim([0 900]); %xlim([0 0.5]);
 ylabel('Ppallet target [Pa]'); xlabel('Schoke/Sj');
 %%
 
 figure();
 scatter(  median( MX(:,:,50)./MX(:,:,6),1,'omitnan'), median(MX(:,:,20)./MX(:,:,19),1,'omitnan') , 'b', 'filled');
 grid on; box on; ylim([0 1]); %xlim([0 0.5]);
 ylabel('Pgroovetonehole/Ppallet');xlabel('Schokevalve/Sj');

 %%
 figure();
 scatter(  median( pi*(0.5*MX(:,:,4)).^2./MX(:,:,6),1,'omitnan'), median(MX(:,:,20)./MX(:,:,19),1,'omitnan') , 'b', 'filled');
 grid on; box on; ylim([0 1]); %xlim([0 0.5]);
 ylabel('Pgroovetonehole/Ppallet');xlabel('Sthonehole/Sj');
%%
 figure();
 scatter(  median(  MX(:,:,5)./MX(:,:,6) ,1,'omitnan'), median(MX(:,:,21)./MX(:,:,19),1,'omitnan') , 'g', 'filled');
 grid on; box on; ylim([0 1]); %xlim([0 0.5]);
ylabel('Pfoot/Ppallet');xlabel('Sin/Sj');
 
 

%%
% ft = fittype( 'Prsv./(1 + ( xx*a + b).^2)','independent','xx','dependent','y');
ft = fittype( 'Prsv./(1 + (a*xx+ b).^2 )','independent','xx','dependent','y');

opts = fitoptions( 'Method', 'NonlinearLeastSquares' , 'TolFun', 1e-12);
opts.Display = 'Off';
opts.Lower      = [Prsv 1e-4 1e-3];
opts.Upper      = [Prsv 5    2]; 
opts.StartPoint = [Prsv 0.5  0.5];
opts.Robust     = 'Bisquare'; % LAR, Off, Bisquare
xData = median( MX(:,:,6)./MX(:,:,5) ,1,'omitnan' );
yData = median( MX(:,:,21)./MX(:,:,19)  , 1,'omitnan');
[fitres, gof] = fit( xData', yData', ft, opts);
figure();
plot(fitres, xData, yData); ylim([0 1]); grid on;xlim([0 1.5]);

%%

Pf_targ    = MX(:,:,21);
Preservoir = 820; % [Pa]

figure(); hold on; grid on; box on;
scatter( sj./sin, Pf_targ, 'b', 'filled'); ylim([0 LM]);

section_sample = [0:1e-5:1.5]; % Sj/Sin
Pfoot_model    = Preservoir./( 1 + (rho/2*section_sample +0.95).^2);

plot(section_sample,Pfoot_model,'--r');
xlabel('Sj/Sin'); ylabel('Pf target [Pa]');



%%


scaledPfootTarg = MX(:,:,21)./MX(:,:,19);

sgr =  PalletDepth*MX(:,:,3);
sin = MX(:,:,5);
sj  = MX(:,:,6);
arearats = (sj./sin).^2 - (sin./sgr).^2;


 figure();
 scatter(arearats  , scaledPfootTarg , 'b', 'filled');
 grid on; ylim([0 1]);

 %%
 zeist = 12*log2(MX(:,:,13)/440);
 figure(); hold on;
 plot(zeist,  log10((sj./sin).^2) , 'bo');
 plot(zeist,  log10((sj./sgr).^2) , 'ro');
 plot(zeist,  log10((sin./sgr).^2) , 'go');



 %%
 alpha = 0.7;
 numer = 1-alpha*xvals.^2;
 denom = 1 + xvals.^2 *(1-alpha);
 figure();
 plot( xvals, numer./denom);




        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Transient analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Assess the goodness of fit for the cases: t80+1PRT, t80+2PRT, t80+3PRT
% figure();
% scatter( 12*log2(MX(:,:,13)/440), MX(:,:,49), 'b', 'filled');
% xlabel('tessitura semitones');
% ylabel('R-squared goodness of fit');
% grid on; box on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Groove 2 Foot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% t20foot - t20groove [groove delay][OK]

figure();
scatter( 12*log2(MX(:,:,13)/440) , 1e3*( MX(:,:,34) - MX(:,:,33)), 'b', 'filled');
xlabel('tessitura [semitones]');
ylabel('t^{20}_{foot} - t^{20}_{groove} [ms]');
box on; grid on; ylim([0 2.2]);

yyaxis right;
scatter( 12*log2(MX(:,:,13)/440), 1e3*MX(:,:,33),'color', "#D95319", 'MarkerFaceColor',	"#D95319");
ylabel('t^{20}{groove} [ms]');
ylim([0 10]);
title('groove-foot delay');


%% PRTfoot vs. PRTgroove [OK]

figure();

scatter( 12*log2(MX(:,:,13)/440), MX(:,:,29)./MX(:,:,28), 'b', 'filled');
ylabel('PRT_{foot}/PRT_{groove} [n.u.]');
grid on; box on;
ylim([0 1.3]);
yyaxis right;
scatter( 12*log2(MX(:,:,13)/440), 1e3*(MX(:,:,29)-MX(:,:,28)), 'color', "#D95319", 'MarkerFaceColor',"#D95319");
ylabel('PRT_{foot} - PRT_{groove} [ms]');
ylim([-1.5 2]);
xlabel('tessitura [semitones]');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% groove-foot analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%>>> tau = Vf/(c_o * S_in)
tau = median(MX(:,:,2)./(co.*MX(:,:,5)), 1, 'omitnan');
tau = 1;
figure();
errorbar(12*log2(F1MEAN/440),...
     median(MX(:,:,34)-MX(:,:,33),1,'omitnan')./tau,...
    std(  (MX(:,:,34)-MX(:,:,33)),1,'omitnan' ) ./ tau,...
    'linestyle','none','color','k','marker', 'o', 'linewidth',1);
ylabel('$t^{20}_{f}-t^{20}_{grv}$','interpreter','latex');
xlabel('Tessitura [semitones]','interpreter','latex');
title('t20 differences');
box on; grid on;
 %ylim([0 2e-3]);

%%
figure();
errorbar( 12*log2(F1MEAN/440),...
    median(MX(:,:,29)./MX(:,:,28),1,'omitnan') ,...
    std( 1*(MX(:,:,29)./MX(:,:,28)),1,'omitnan' ),...
    'linestyle','none','color','k','marker', 'o', 'linewidth',1);

ylabel('$PRT_{f}/PRT_{grv}$','interpreter','latex'); 
xlabel('Tessitura [semitones]','interpreter','latex');
title('PRT ratios');
grid on; box on; ylim([0 1.3]);



%% PRTfoot vs Vf

figure();
errorbar( median(MX(:,:,2),1,'omitnan'),...
    median(MX(:,:,29),1,'omitnan'),...
    std(MX(:,:,29)),...
    'linestyle','none','color','k','marker', 'o', 'linewidth',1);





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transient spectrum study
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% max(a2/a1) during transient [not bad]

% GVC: I have used butter(4), filtfilt() of envel_first and envel_second
% and then search max within t20_f and t80_f+50*PRT.

mask = [1:22]; %mask(7)=[]; % >>> If you want to remove a single dataset

figure();
% scatter( 12*log2(MX(:,mask,13)/440), MX(:,mask,48), 'b', 'filled');
errorbar( 12*log2(median(MX(:,mask,13)/440,1,'omitnan')), median(MX(:,mask,48),1,'omitnan') , std(MX(:,mask,48),1,'omitnan') );
%%
grid on; box on;
xlabel('semitones wrt 440');ylabel('max( a_2(t)/a_1(t) smooth) during trans. (lin) [n.u.]');
hold on;
plot( 12*log2(F1MEAN(mask)/440),median(MX(:,mask, 48),1,'omitnan') ,  '-ro');



figure();
boxplot(MX(:,mask,48), 'Positions', 12*log2(F1MEAN(mask)/440), 'Widths',1);
box on; grid on;
xlabel({'xticks sample pipe','spacing semitones wrt 440'});
ylabel('max( a_2(t)/a_1(t) smooth) during trans. [n.u.]'); 
ylim([0 50]);



%% Transient spectra differences integrated in area between a_2(t) and a_1 for t=[t20_p, t80_p] [not bad]

figure();
scatter( 12*log2(MX(:,:,13)/440),    MX(:,:,38).*MX(:,:,13)   ,'b' , 'filled');
hold on;
plot( 12*log2(F1MEAN/440), median( MX(:,:,38).*MX(:,:,13), 1, 'omitnan')  , '-ro');
ylabel('T1 * trapez(dt, a_2(t) - a_1(t) )  ,  t=[t^{20}_{p}, t^{80}_{p}] [n.u.]'); xlabel('tessitura semitones');
xlabel('tessitura');
grid on; box on;
title('Integrated area comprised between a_2(t) and a_1(t) in the interval (t^{20}_p, t^{80}_p) norm. by oscill. period');
ylim([-5 40]);

%% Transient spectra differences integrated in area between a_3(t) and a_1 for t=[t20_p, t80_fp [not bad]

figure();
scatter( 12*log2(MX(:,:,13)/440),    MX(:,:,39).*MX(:,:,13)   ,'b' , 'filled');
hold on;
plot( 12*log2(F1MEAN/440), median( MX(:,:,39).*MX(:,:,13), 1, 'omitnan')  , '-ro');

ylabel('T1 * trapez(dt, a_3(t) - a_1(t) ), t=[t^{20}_p, t^{80}_p)] [n.u.]'); xlabel('tessitura [semitones]');

grid on; box on;
title('Integrated area comprised between a_3(t) and a_1(t) in the interval (t^{20}_p, t^{80}_p) norm. by oscill. period');
ylim([-5 40]);


%% Fig6, CFA Ernoult2016:: Nondim. of a2 versus \theta at t= t(a_{2max}) [meeeh...]

mask = [1:22]; mask([1,2])=[]; % <<<<< If you want to filter out some sample pipes

Uj_at_a2max =  sqrt(MX(:,mask,46)*2/1.2);
theta_a2max = Uj_at_a2max./(MX(:,mask,9).*MX(:,mask,13));
a2_nondim   = MX(:,mask,47)./(rho*co*Uj_at_a2max);

figure();
% scatter( 12*log2(MX(:,:,13)/440),a2_nondim);
scatter(theta_a2max, a2_nondim,'filled');
xlabel('\theta at t=a_{2,f_2max} [n.u.]');
ylabel('a_{2,max} / \rho c_o U_j  @ t(a_{2,max}) [n.u.]');
grid on; box on;
xlim([6 11]);ylim([0 2.2e-3]);



%% a2max over a2 target [OK]
figure();
scatter( 12*log2(MX(:,:,13)/440), MX(:,:,45), 'b', 'filled');
xlabel('Tessitura [semitones]');
ylabel('a_{2,max} / a^{targ}_{2}  [n.u.]');
grid on; box on;ylim([0 5]);

%% a2max/a1(@a2max) (mouth rad) [OK]
figure();
scatter(12*log2(MX(:,:,13)/440), max(0,MX(:,:,43)), 'b','filled' );
ylabel('a_{2,max}/a_{1}^{simult.} [n.u.]'); 
box on; grid on; ylim([0 25]);

figure();
funh = boxplot( (abs(MX(:,:,43))) , 'Positions',12*log2(F1MEAN/440), 'Widths',1);
ylabel('a_{2,max} / a_{1}^{simult.} [n.u.]');
ylim([0 25]); grid on;
xlabel({'Sample pipe xticks','12\times  log_2(f_1/440) spacing'},'interpreter','latex');




%% a2max/a1target (mouth rad) [OK]

figure();
scatter(12*log2(MX(:,:,13)/440), (MX(:,:,44)), 'b','filled' );
xlabel('$12\times log_2(f_1/f_{ref})$','interpreter','latex');
ylabel('$a_{2,max}/a_1^{targ]}$ (lin)','interpreter','latex'); 
box on; grid on; ylim([0 2.2]);

figure();boxplot( (MX(:,:,44)), 'Positions',12*log2(F1MEAN/440), 'Widths',1);
xlabel({'sample pipe xticks','12\times log_2(f1/440) spacing'},'interpreter','latex');
ylabel('$a_{2,max}/a_{1}^{targ}$ (lin)','interpreter','latex');
 ylim([0 2.2]); grid on;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ????????
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pmouth as per massage to equations 9-10-11 and alpha_vc = 1;


Pm =  MX(:,:,21) - (MX(:,:,5)./MX(:,:,6)).^2.*(MX(:,:,20) - MX(:,:,21));
Pm =  MX(:,:,21) + (MX(:,:,5)./MX(:,:,6)).^2.*(MX(:,:,21) - MX(:,:,20));

figure();

scatter(  12*log2(MX(:,:,13)/440), log10(Pm) , 'filled'); box on; grid on;

xlabel('$12\times  log_2(f_1/440)$','interpreter','latex'); box on; grid on;
ylabel('$P_{m}$ [Pa] ($log_{10}$)','interpreter','latex');
title('Expected mouth-rad pressure as per eqs. 9-10-11 of model and Qj=Qin, no alpha_vc');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delays (t20): foot and mouth rad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% t20grv - t20 foot [ACOUSTIC DELAY][OK] 06/12/2024[[[[[[06/12/2024 [to keep]]]]]]]

figure();
scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan'), median(1e3*abs(MX(:,:,33)-MX(:,:,34) ),1,'omitnan') , 'b', 'filled');
ylabel('$t^{20}_{grv}$ - t$^{20}_{foot}$ (lin)[ms]','interpreter','latex');
xlabel('$12log_2(f_1/440)$','interpreter','latex'); box on; grid on;

ylim([0 2]);

%% t80grv - t80 foot [ACOUSTIC DELAY][OK] [[[[[[06/12/2024 [to keep]]]]]]]

figure();
scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan'), ...
    median(   MX(:,:,33)+MX(:,:,28)-MX(:,:,29)-MX(:,:,34) ,1,'omitnan') ,...
    'b', 'filled');

ylabel('$t^{20}_{grv}$ - t$^{20}_{foot}$ (lin)[ms]','interpreter','latex');
xlabel('$12log_2(f_1/440)$','interpreter','latex'); box on; grid on;

%% t80grv - t80 foot [ACOUSTIC DELAY][OK] 06/12/2024 [to keep][[[[[[06/12/2024 [to keep]]]]]]]

figure();
scatter( median( MX(:,:,2)./MX(:,:,11),1,'omitnan'), ...
    median(  1e3* abs(MX(:,:,33)+MX(:,:,28)-MX(:,:,29)-MX(:,:,34)) ,1,'omitnan') ,...
    'b', 'filled');







%%
figure();
scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan'), median(1e3*abs(MX(:,:,28) ),1,'omitnan') , 'b', 'filled');
title(' PRT groove');
ylim([0 4]);
%%
figure();
scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan'), median(1e3*abs(MX(:,:,29) ),1,'omitnan')./median(1e3*abs(MX(:,:,28) ),1,'omitnan') , 'b', 'filled');
title('PRT foot');


ylim([0 1.3]);



%% t20rad - t20 foot [ACOUSTIC DELAY][OK]

figure();
scatter( 12*log2(MX(:,:,13)/440), 1e3*(MX(:,:,35)-MX(:,:,34) ) , 'b', 'filled');
ylabel('$t^{20}_{rad}$ - t$^{20}_{foot}$ (lin)[ms]','interpreter','latex');
xlabel('$12log_2(f_1/440)$','interpreter','latex'); box on; grid on;
title('Acoustic delay','interpreter','latex');ylim([0 80])

%% (t20rad-t20foot)/T1 [ACOUSTIC DELAY][OK]

figure();
scatter( 12*log2(MX(:,:,13)/440), (MX(:,:,35)-MX(:,:,34) ).*MX(:,:,13), 'b', 'filled');
ylabel('$( t^{20}_{rad} - t^{20}_{foot} ) /T_1 $','interpreter','latex');%, 'Rotation',0);
xlabel('$12log_2(f_1/440)$','interpreter','latex'); box on; grid on;
title('Acoustic delay','interpreter','latex');%ylim([0 ])

%% t20_foot [HYDRODYNAMIC DELAY][OK]

figure();
scatter( 12*log2(MX(:,:,13)/440),1e3*(MX(:,:,34) ), 'b', 'filled');
ylabel('t$^{20}_{foot}$ (lin)[ms]','interpreter','latex');
xlabel('$12log_2(f_1/440)$','interpreter','latex'); box on; grid on;
title('Hydrodynamic delay','interpreter','latex');ylim([0 10])

%% t20_foot normalized by tau_f = Vf/(c_o * S_in) [HYDRODYNAMIC DELAY]

figure();
scatter( 12*log2(MX(:,:,13)/440),( co * MX(:,:,34).*MX(:,:,5)./MX(:,:,2) ), 'b', 'filled');
ylabel('$ t^{20}_f / \tau_f $ with $\tau_f = V_f / c_o  \mathcal{S}_{in}$','interpreter','latex');
xlabel('$12\times log_2(f_1/440)$','interpreter','latex'); box on; grid on;
title('Hydrodynamic delay','interpreter','latex');
ylim([0 1.3]);

%% t20_foot normalized by tau_f = Vf/(u_in * S_in) [HYDRODYNAMIC DELAY]

tau_f = MX(:,:,2)./(MX(:,:,5).*sqrt( 2/rho*(MX(:,:,20)-MX(:,:,21))) );

figure();
scatter( 12*log2(MX(:,:,13)/440), MX(:,:,34)./tau_f , 'b', 'filled');
ylabel('$ t^{20}_f / \tau_f $ with $\tau_f = V_f / u_{in}  \mathcal{S}_{in}$','interpreter','latex');
xlabel('$12\times log_2(f_1/440)$','interpreter','latex'); box on; grid on;
title('Hydrodynamic delay','interpreter','latex');
ylim([0 0.2]);

%% t20_foot normalized by tau_f = Vf/(u_j * S_j) [HYDRODYNAMIC DELAY][interesting...]

tau_f = MX(:,:,2)./(MX(:,:,6).*sqrt( 2/rho*(MX(:,:,21)-0*MX(:,:,21))) );

figure();
scatter( 12*log2(MX(:,:,13)/440), MX(:,:,34)./tau_f , 'b', 'filled');
xlabel('Tessitura [semitones]');
ylabel('$ t^{20}_f / \tau_f $ with $\tau_f = V_f / u_{j}  \mathcal{S}_{j}$','interpreter','latex');
grid on; box on;
ylim([0 0.055]);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRT's 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PRTfoot/T1 (w.r.t tessitura) [T1 for normalization??? it's hydrodynamic]
figure();

tau_f = MX(:,:,2)./(MX(:,:,5).*sqrt( 2/rho*(MX(:,:,20)-MX(:,:,21))) );

% scatter(12*log2(MX(:,:,13)/440), (MX(:,:,29).*MX(:,:,13)), 'b','filled'); ylim([0 2]);
scatter(12*log2(MX(:,:,13)/440), MX(:,:,29)./tau_f, 'b','filled' ); ylim([0 0.028]);
xlabel('$12\times log_2(f_1/f_{ref})$','interpreter','latex');
ylabel('$\frac{PRT_{foot}}{T_1}$','interpreter','latex','Rotation',0);
box on; grid on;

%% PRTrad/T1 (w.r.t tessitura) [OK]
figure();
errorbar(    median(12*log2(MX(:,:,13)/440),1,'omitnan'), ...
    median(MX(:,:,30).*MX(:,:,13), 1, 'omitnan'),...
    std(MX(:,:,30).*MX(:,:,13), 'omitnan'), ...
    'linestyle','none', 'marker', 'o', 'color', 'k', 'linewidth', 1); 



xlabel('$12 \times log_2(f_1/f_{ref})$','interpreter','latex');
ylabel('$ PRT_{rad} /T_1$','interpreter','latex');
box on; grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Study of Beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% beta vs t20foot [meeh...]

figure();

scatter( MX(:,:,26), 1e3*MX(:,:,34), 'b', 'filled');
grid on; box on; ylim([0 11]);
xlabel('beta [s^{-1}]'); ylabel('t20foot [ms]')

%% beta vs PRTfoot [OK]
figure();
scatter( 1./MX(:,:,26), 1e3*MX(:,:,29)  , 'b','filled');
xlabel('$\beta^{-1} [s]$', 'interpreter','latex');
ylabel('PRT$_{foot}$ [ms]', 'Interpreter','latex');
box on; grid on;
xlim([0 2.2e-3]);
ylim([0 4.5]);

%% nu w.r.t. tessitura [ok]

figure();
scatter( 12*log2(MX(:,:,13)/440) , MX(:,:,27), 'b', 'filled');
xlabel('tessitura [semitones]');
ylabel('nu [n.u.]');
grid on; box on; 
hold on;
plot( 12*log2(F1MEAN/440) , median(MX(:,:,27), 1, 'omitnan') , '-ro');

%% BETA*T1 versus NU [OK]
figure();
scatter( (MX(:,:,26)./MX(:,:,13)), (MX(:,:,27)), 'b','filled');
xlabel('beta*T1 [n.u.]'); ylabel('nu [n.u.]');
grid on; box on; axis equal;

%% beta vs exp(nu) [maybe?]
figure();
mask = [1:10,13:22];
mask = 1:22;

scatter( (MX(:,mask,26)), (MX(:,mask,27)), 'b', 'filled');
xlabel('$\beta [s^{-1}]$', 'interpreter','latex');
ylabel('$\nu [n.u.]$', 'Interpreter','latex');
grid on; box on;

%% nu vs beta^(0.56) [meeh...]
figure();
mask = [1:10,13:22];
mask = 1:22;
% scatter( MX(:,mask,26), (MX(:,mask,27)), 'b', 'filled');
% scatter( log10(MX(:,mask,26))/log10(20), log(MX(:,mask,27)), 'b', 'filled');
scatter( MX(:,mask,27), MX(:,mask,26).^(0.56), 'b', 'filled');
ylabel('$\beta^{0.56}$ [units?]', 'interpreter','latex');
xlabel('$\nu$', 'Interpreter','latex'); grid on; box on;
ylim([0 90]);


%%  Beta x T1 vs freq [OK]

mask = [1:22]; % low: sample pipes 2 and 3; high, sample pipes 11 and 12

figure();
% scatter( 12*log2(MX(:,mask,13)/440) , log10(MX(:,mask,26)./MX(:,mask,13)) ,'b','filled');ylabel('beta x T_1 (log10)');xlabel('12log_2(f_1/440)'); box on; grid on; 
scatter( 12*log2(MX(:,mask,13)/440) , (MX(:,mask,26)./MX(:,mask,13)) ,'b','filled');
ylabel('$\beta \times  T_1$ [n.u.]','interpreter','latex');
xlabel('$12log_2(f_1/f_{ref})$','interpreter','latex'); 
box on; grid on; 
ylim([0 11]);
hold on;
plot( 12*log2(F1MEAN/440), median(MX(:,mask,26)./MX(:,mask,13), 1, 'omitnan')   , '-or');

%% Beta x T1 vs freq BOX plot [OK]

figure();
boxplot(MX(:,:,26)./MX(:,:,13), 'Positions',12*log2(F1MEAN/440), 'Widths',1);
xlabel({'Sample pipe xticks','tessitura fspacing'},'interpreter','latex');
ylabel('$\beta \times T_1$','interprete','latex');
ylim([0 11]); grid on;

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

%%
figure();
scatter( 12*log2( MX(:,:,13)/440 ), log10(MX(:,:,17)) ,'filled');


%% PRTgrv * f1 vs NU [NO]
figure();
scatter( log10(MX(:,:,30).*MX(:,:,13)), log10(MX(:,:,27)), 'b', 'filled');
%% PRTgrv * f1 vs NU
figure();
scatter( log10(MX(:,:,30).*MX(:,:,13)), log10(MX(:,:,27)), 'b', 'filled');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Old vena-contracta investigations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% vena-contracta factor groove [OK][But what does it mean?]

FSZ = 14;

factor = (MX(:,:,5)./(MX(:,:,3)*PalletDepth)).*sqrt( (MX(:,:,20) - MX(:,:,21) )./(MX(:,:,19)-MX(:,:,20)));

figure();
scatter( MX(:,:,5)*1e6, factor, 'b', 'filled');box on;
xlabel('$S^{geo}_{in} \ [mm^2]$','interpreter','latex','fontsize',FSZ);
ylabel('$\Gamma_{grv}$','interpreter','latex','fontsize',FSZ);


%% #################################################################### 
% Ratio of effective areas as a function of geometric areas ratio
% #################################################################### 

FSZ = 14;
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





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% PREPARE PCA()  %%%%%%%%%%%%%%%%%%%%

MXresh = nan*ones(size(MX,1)*size(MX,2),size(MX,3));
for idx = 1 : size(MX,3)
   MXresh(:,idx) = reshape( MX(:,:,idx),numel(MX(:,:,idx)),1 ) ;  
end

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
scatter( freqlogax, log10(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
hold on;
yyaxis right;
scatter( freqlogax, log10(MX(:,:,27)), 'r', 'filled');box on;
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
scatter( freqlogax, log10(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
hold on;
yyaxis right;
scatter( freqlogax, log10(MX(:,:,27)), 'r', 'filled');box on;
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
scatter( log10(MX(:,:,38)), log10(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
xlabel('I21');

figure();
scatter( log10(MX(:,:,38)), log10(MX(:,:,27)), 'r', 'filled');box on;
ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
xlabel('I21');

%% 
% #################################################################### 
% BETA and NU with PTARG GROOVE
% #################################################################### 


figure(27); clf;
FSZ = 15;
scatter( log10(MX(:,:,20)), log10(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
hold on;
yyaxis right;
scatter( log10(MX(:,:,20)), log10(MX(:,:,27)), 'r', 'filled');box on;
% ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
% legend('$\beta$','$\nu$','interpreter','latex','fontsize',FSZ);
yyaxis right;
ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
xlabel('$log10 PTARG groove$','interpreter','latex','fontsize',FSZ);



%% 
% #################################################################### 
% BETA and NU with PTARG FOOT
% #################################################################### 


figure(28); clf;
FSZ = 15;
scatter( log10(MX(:,:,21)), log10(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
hold on;
yyaxis right;
scatter( log10(MX(:,:,21)), log10(MX(:,:,27)), 'r', 'filled');box on;
% ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
% legend('$\beta$','$\nu$','interpreter','latex','fontsize',FSZ);
yyaxis right;
ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
xlabel('$log10 PTARG foot$','interpreter','latex','fontsize',FSZ);


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
scatter( freqlogax, log10(MX(:,:,1)), 'b', 'filled');box on;
ylabel('$Lp$','interpreter','latex','fontsize',FSZ);

geoh(2) = subplot(3,4,2);
scatter( freqlogax, log10(MX(:,:,2)), 'b', 'filled');box on;
ylabel('$Vf$','interpreter','latex','fontsize',FSZ);

geoh(3) = subplot(3,4,3);
scatter( freqlogax, log10(MX(:,:,3)), 'b', 'filled');box on;
ylabel('$PWdth$','interpreter','latex','fontsize',FSZ);

geoh(4) = subplot(3,4,4);
scatter( freqlogax, log10(MX(:,:,4)), 'b', 'filled');box on;
ylabel('$TnHD$','interpreter','latex','fontsize',FSZ);

geoh(5) = subplot(3,4,5);
scatter( freqlogax, log10(MX(:,:,5)), 'b', 'filled');box on;
ylabel('$Sin$','interpreter','latex','fontsize',FSZ);

geoh(6) = subplot(3,4,6);
scatter( freqlogax, log10(MX(:,:,6)), 'b', 'filled');box on;
ylabel('$Sjet$','interpreter','latex','fontsize',FSZ);

geoh(7) = subplot(3,4,7);
scatter( freqlogax, log10(MX(:,:,7)), 'b', 'filled');box on;
ylabel('$h$','interpreter','latex','fontsize',FSZ);

geoh(8) = subplot(3,4,8);
scatter( freqlogax, log10(MX(:,:,8)), 'b', 'filled');box on;
ylabel('$H$','interpreter','latex','fontsize',FSZ);

geoh(9) = subplot(3,4,9);
scatter( freqlogax, log10(MX(:,:,9)), 'b', 'filled');box on;
ylabel('$Wm$','interpreter','latex','fontsize',FSZ);

geoh(10) = subplot(3,4,10);
scatter( freqlogax, log10(MX(:,:,10)), 'b', 'filled');box on;
ylabel('$Dp$','interpreter','latex','fontsize',FSZ);

% geoh(11) = subplot(3,4,11);
% scatter( freqlogax,log10( MX(:,:,11)), 'b', 'filled');box on;
% ylabel('$Vgrv$','interpreter','latex','fontsize',FSZ);

geoh(12) = subplot(3,4,12);
scatter( freqlogax, log10(MX(:,:,12)), 'b', 'filled');box on;
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
scatter( freqlogax, log10(MX(:,:,26)), 'b', 'filled');box on;
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
scatter( freqlogax, log10(MX(:,:,38)), 'b', 'filled');box on;ylim([-9 -1]);
ylabel('$I1 \ (log10)$','interpreter','latex','fontsize',FSZ);

f2h(19)=subplot(8,3,13);
scatter( freqlogax, MX(:,:,24), 'b', 'filled');box on;
ylabel('$P^{o}_{(ft-grv)}$','interpreter','latex','fontsize',FSZ);

f2h(20)=subplot(8,3,17);
scatter( freqlogax, MX(:,:,32), 'b', 'filled');box on;
ylabel('$PRT_{(rad/ft)}$','interpreter','latex','fontsize',FSZ);

f2h(21)=subplot(8,3,[21,24]);
scatter( freqlogax, log10(MX(:,:,39)), 'b', 'filled');box on;ylim([-9 -1]);
ylabel('$I2 \ (log10)$','interpreter','latex','fontsize',FSZ);

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

% #################################################################### 
% Q-factor vs BETA
% #################################################################### 

figure(12);clf; hold on;
% Q-factor vs Beta
try
for idx = 1 : length(QFAC1)
    scatter(QFAC1(idx), log10(MX(:, idx, 5)), 'b','filled');
end
end
xlabel('$Q-$factor [s$^{-1}$]','interpreter','latex','fontsize',14);
ylabel('$\beta$ fitted (log10) [s$^{-1}$]','interpreter','latex','fontsize',14);
hold on

BM = MX(:, :, 5);
BM = mean(BM,1, 'omitnan');
p = polyfit(QFAC1(find(BM)), log10(BM), 1);
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
scatter(  log10(Vf_mask(idx)), log10(NU(:,idx)), 'b','filled'); % YES good
% scatter( log10( ratios(idx) ), log10((NU(:,idx))), 'b','filled'); 
% scatter( log10( Inletmask(idx)./palletarea_mask(idx) ), log10((NU(:,idx))),'b','filled'); 
% scatter( log10( Sjetmask(idx)./palletarea_mask(idx) ), log10((NU(:,idx))),'b','filled'); 
% scatter( log10( Inletmask(idx) ), log10((NU(:,idx))),'b','filled'); 

end
end


hold on;

NM = mean(NU,1, 'omitnan');
p = polyfit(log10(ratios), log10(NM), 1);
nulin = polyval(p, log10(ratios) );
% plot(log10(ratios), (nulin), '--r')


FSZ = 14;
xlabel('$S_{jet}/S_{inlet}$ (log10)','interpreter','latex','fontsize',FSZ);
ylabel('$\nu$ (log10)','interpreter','latex','fontsize',FSZ);
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
