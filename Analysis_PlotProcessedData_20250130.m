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

% [28]:PRT20grv               [29]:PRT20foot              [30]:PRT20pipe
% [31]:PRT20foot/grv          [32]:PRT20rad/foot  
% [33]:t20 groove           [34]:t20 foot             [35]:t20 rad
% [36]:t20 foot-groove      [37]:t20 rad-foot
% [38]:Area1                [39]:Area2
% [40]:Sin/Spall            [41]:Sjet/Sin             [42]: KeyVel

% [43]: A2max_over_A1simult [44]: A2max_over_A1target [45]: a2max_over_a2target
% [46]: DeltaP(foot-mouth)_at_a2max
% [47]: a2max_vec        [48]: max_a2_over_a1 (after smooth)\in(t20_f,t80_f+50PRT)
% [49]: gofr2 (r-squared goodness of logistic fit)
% [50] Lateral stroke area of pallet valve at max opening (smaller than Slot Area = PWidth*Length_win_slot =MX(:,:,3)*0.129;


% [51]: t10groove [52]: PRT10groove 
% [55]: t10foot   [56]: PRT10foot 
% [57] t10pipe   [58]: PRT10pipe


rho = 1.2;
co  = 340;

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %       Geometry
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% All geometry vars w.r.t. f1/440
FSZ = 17;
figure(); hold on;

for idx = 1 :12
    tmp_data = log10( median(MX(:,:,idx)) );
    val_absc = interp1(freqlogax, tmp_data, 0);
    plot(freqlogax, tmp_data-val_absc ,'-') ;
    
end
grid on; box on;



namevarsall = {'Lp','Vf','PWdth','TNHD','INLET','SJET','h','H','Wm','Dp',...
    'Vgrv','Qfact','f1','theta','Qpall2grv','Qgr2ft','Qjet','Rmplssg',...
    'P{o} pall','P{o}grv','P{o}foot','P{o}rad',...
    'P{o} grv-pall','P{o}foot-grv','P{o}rad-foot',...
    'Beta','nu',...
    'PRTgrv','PRTfoot','PRTrad',...
    'PRTfoot/grv','PRTrad/foot',...
    't{20}grv','t{20}foot','t{20}rad',...
    't20(foot-grv)','t20(rad-foot)',...
    'Area1','Area2','Sin/Spal','Sjet/Sin','MaxKeyVel',...
    'a2maxOverA1simult','a2maxOvera1Targ','a2maxOvera2Targ','DeltaPfoot2mouthAta2max','a2maxVec','a2maxOverA1smooth',...
    'gofr2','LateralChkSec'};
legend(namevarsall);

figure(); 
for idx=1:12

    ff = polyfit( 12*log2(F1MEAN/440), log10(median(MX(:,:,idx),1,'omitnan')), 1);
    plot(12*log2(F1MEAN/440), log10(median(MX(:,:,idx),1,'omitnan')) ,'-kd');
    hold on;
    plot( 12*log2(F1MEAN/440), polyval(ff,12*log2(F1MEAN/440)),'--r');
    grid on; box on;
    hold off;
    title(sprintf('Variable: %s. Slope %2.5f, offset %2.5f ',namevarsall{idx}, ff(1),ff(2)));
    drawnow;
    pause();
end




%% Sj / Sin [OK]

figure();
scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan'), median((MX(:,:,6)./MX(:,:,5)),1,'omitnan'), 'kd', 'filled');
ylabel('$\mathcal{S}_j / \mathcal{S}_{in} \ (log_{10})$ [n.u.]', 'interpreter','latex');
xlabel('$12log_2(f_1/440)$','interpreter','latex'); box on; grid on; ylim([0 1.5]);




%% Sj / Sin */ Vf^0.33 [interactions]

figure();
scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan'), median(   (MX(:,:,6)./MX(:,:,5))./MX(:,:,2).^(1/3)   ,1,'omitnan'), 'kd', 'filled');
ylabel('$\mathcal{S}_j / \mathcal{S}_{in} /Vf^{(1/3)}\ (log_{10})$ [n.u.]', 'interpreter','latex');
xlabel('$12log_2(f_1/440)$','interpreter','latex'); box on; grid on;ylim([0 20]);


%% S_palletvalve, S_slot, S_tonehole, S_in, S_j 


figure(); hold on; grid on; box on;

scatter(freqlogax, median(1e6*MX(:,:,50),1,'omitnan'), 'b', 'filled');       % mean = 8.4 cm^2 % pallet valve stroke
scatter(freqlogax, median(1e6*MX(:,:,3)*0.1298,1,'omitnan'), 'r', 'filled'); % mean = 16.3 cm^2 % pallet slot window
scatter(freqlogax, median(1e6*pi*(0.5*MX(:,:,4)).^2,1,'omitnan'),'k','filled' ); % tone hole area
scatter(freqlogax, median(1e6*MX(:,:,5),1,'omitnan'), 'g', 'filled'); % foot inlet
scatter(freqlogax, median(1e6*MX(:,:,6),1,'omitnan'), 'm', 'filled'); % foot outlet
ylabel('mm^2'); xlabel('tessitura [semitones]');
ax=gca; ax.YScale = 'log';





                            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %       Steady-State analysis
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Estimation of groove flow velocity


Sgroove = 0.05*MX(:,:,3); % No vena contracta effect, but boundary layer that we'll neglect for now
ug = 1.5*MX(:,:,6)./Sgroove.*sqrt(2/rho.*MX(:,:,21));

figure();
scatter( 12*log2(F1MEAN/440), median(ug,1,'omitnan'), 'kd','filled');

        
%% Pressure drops versus cross sections [not yet there...]
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
box on; grid on;

%% Foot flow conservation and Gamma function [OK][Keep, 2025/01/30, plot1/2 venacontracta]
      
Qin = 1.0*MX(:,:,5) .* sqrt( 2/rho*( MX(:,:,19) - MX(:,:,21) ) );
Qj  = 1.0*MX(:,:,6).*sqrt( 2/rho * ( MX(:,:,21) -      0     ) );
 
FitFlowConserv = polyfit( median(MX(:,:,6)./MX(:,:,5),1,'omitnan') , median(Qj./Qin,1,'omitnan') ,1);

figure();
errorbar(  median(MX(:,:,6)./MX(:,:,5),1,'omitnan') , median( Qj./Qin,1,'omitnan'), std( abs(Qj)./Qin,1, 'omitnan'),'k', 'marker','d','linestyle','none');
% errorbar(  median(12*log2(MX(:,:,13)/440),1,'omitnan') , median( abs(Qj)./Qin,1,'omitnan'), std( abs(Qj)./Qin,1, 'omitnan'),'k', 'marker','x','linestyle','none');
grid on; ylim([0 1]); hold on;
xlabel('$ \mathcal{S}_j / \mathcal{S}_{in}$','interpreter','latex');
ylabel('$\Gamma (\mathcal{S})$','interpreter','latex');

querypoints = [0:1e-2:1.5];
plot(querypoints, polyval(FitFlowConserv, querypoints), '--k');

text(0.8,0.55,sprintf('$y=%1.3f x + %1.3f$',FitFlowConserv(1),FitFlowConserv(2) ), 'interpreter','latex');


%% Theta
figure();

scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan'), median(MX(:,:,14),1,'omitnan'), 'dk', 'filled');
grid on; box on;
xlabel('$12 \times log_2 (f_1 /  440 Hz) $', 'interpreter','latex');
ylabel('$\theta = u_j/f_1 W_m$','interpreter','latex'); ylim([0 12]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target pressures study
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Ptarg differences all wrt to reservoir, in SS grv and rsv have the same pressure [OK,final][keep][2025/01/30]
figure();hold on;
    plot( 12*log2(F1MEAN/440), median(abs(MX(:,:,19)-MX(:,:,20))./MX(:,:,19) , 1, 'omitnan'),...       
        'k','linestyle','none','marker','o');
    plot( 12*log2(F1MEAN/440), median(abs(MX(:,:,20)-MX(:,:,21))./MX(:,:,19) , 1, 'omitnan'),...
        'k','linestyle','none','marker','^');
    plot( 12*log2(F1MEAN/440), median(abs(MX(:,:,21)-0)./MX(:,:,19) , 1, 'omitnan'),...
        'k','linestyle','none','marker','x');

xlim([-20 25]);ylim([0 1]); grid on; box on;
xlabel('$12\times log_2 (f_1/440 \ Hz)$','interpreter','latex');
legend('$\widetilde{\Delta P}_{pall2grv} = |P^{\oplus}_{gr}-P^{\oplus}_{pall}|/P^{\oplus}_{pall}$',...
    '$\widetilde{\Delta P}_{grv2ft}\  = |P^{\oplus}_{ft}-P^{\oplus}_{grv}|/P^{\oplus}_{pall}$',...
    '$\widetilde{\Delta P}_{ft2m}\ \   = |\langle P^{\oplus}_{m}\rangle-P^{\oplus}_{ft}|/P^{\oplus}_{pall}$','interpreter','latex');

%% Pfoot target versus Sj/Sin [ok][keep][2025/01/30][plot2/2, venacontracta]

 % [Can we predict PfootTarg just by Sj and Sin?]
 
 Prsv   = 820.83;
 Prsv   = 1;
 S_axes = linspace(0,1.5, 1e3);
 
 % [1] DATA measured: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Sgeom_ratio_measured = median(MX(:,:,6)./MX(:,:,5),1,'omitnan');
 Qj_over_Qin    = median(  (MX(:,:,6).*sqrt(2/rho*(MX(:,:,21)-0)))  ./  (MX(:,:,5).*sqrt(2/rho*(MX(:,:,20)-MX(:,:,21)))) ,...
           1,'omitnan');
 
 
FitFlowConserv = polyfit( Sgeom_ratio_measured , Qj_over_Qin ,1);
  
A = FitFlowConserv(1); 
B = FitFlowConserv(2);
Gamma    = (A * S_axes + B); 
GammaInv = 1./Gamma;
 
Sgeom_measured = median(MX(:,:,6)./MX(:,:,5),1,'omitnan');
Seff_measured  = 0.5477*Sgeom_measured + 0.1398;
 
% gamma_measured = median(MX(:,:,6)./MX(:,:,5).*sqrt(MX(:,:,21)./abs(MX(:,:,21)-MX(:,:,19))),1,'omitnan');
gamma_measured  = Qj_over_Qin;
 
 
% [2] MODEL:  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
Pfoot_model = Prsv./(1 + (GammaInv.*S_axes).^2);
x = linspace(0,2.0,1e2);

figure(); hold on;
if 0           % NO:
 plot( 1./gamma_measured .* Sgeom_ratio_measured, median(MX(:,:,21)./MX(:,:,19),1,'omitnan') , 'xk');
 y = 1./(1+x.^2);
end
if 0            % Poor:
    plot( 1./(A * Sgeom_ratio_measured + B).*Sgeom_ratio_measured , median(MX(:,:,21)./MX(:,:,19),1,'omitnan') , 'dk');
    y = 1./(1+x.^2); % <<< ??????
end
if 1            % Best:
    plot( Sgeom_ratio_measured , median(MX(:,:,21)./MX(:,:,19),1,'omitnan') , 'dk','markerfacecolor','k');
 y = 1./(1 + (x./(A*x+B)).^2);
 xlim([0 1.5]);
end
 
 
 plot(x, y, '--k');
 ylim([0 1]);
 
 grid on; box on; 
 ylabel('$P_f^{\oplus}/P_{pall}^{\oplus}$ [n.u.]', 'interpreter','latex');
 xlabel('$\frac{\mathcal{S}_j}{ \mathcal{S}_{in}}$ [n.u.]', 'interpreter','latex');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uj/uin (no vena contracta) [OK]

figure();
scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan') , median(MX(:,:,18),1,'omitnan'), 'dk', 'filled');
grid on; box on;
ylim([0 1]);
xlabel('Tessitura [semitones]'); title('Remplissage');
 ylabel('$u_j / u_{in} = \sqrt{P_f^o} / \sqrt{ P_{gr}^o - P_f^o} $','interpreter','latex');

 % Additional hypothesis on why we don't see the flow consrevation:
 % the foot is a big volume with a small channel at the entrance and output
 % so we have a jet at the entrance dissipated by turbulence
 % and a second at the output....? so?


%%
ft = fittype( 'Prsv./(1 + (a*xx + b).^2 )','independent','xx','dependent','y');

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

%% ????? Pmouth as per massage to equations 9-10-11 and alpha_vc = 1;

Pm =  MX(:,:,21) - (MX(:,:,5)./MX(:,:,6)).^2.*(MX(:,:,20) - MX(:,:,21));
Pm =  MX(:,:,21) + (MX(:,:,5)./MX(:,:,6)).^2.*(MX(:,:,21) - MX(:,:,20));

figure();
scatter(  12*log2(MX(:,:,13)/440), log10(Pm) , 'filled'); box on; grid on;
xlabel('$12\times  log_2(f_1/440)$','interpreter','latex'); box on; grid on;
ylabel('$P_{m}$ [Pa] ($log_{10}$)','interpreter','latex');
title('Expected mouth-rad pressure as per eqs. 9-10-11 of model and Qj=Qin, no alpha_vc');



                        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %       Transient analysis
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Characteristic times ============

Vgrv = median(MX(:,:,11),1,'omitnan');
Vf   = median(MX(:,:,2),1, 'omitnan');

Spall_Slot    = median(MX(:,:,3),1,'omitnan')*0.1298; % Perforated rectangles on the plate, with the same width as the groove channel
Spall_Lateral = PalletValveStrokeArea(maskpipes); % At maximum aperture of valve, adding areas of a rectangle and two triangles
Sin = 1*median(MX(:,:,5),1,'omitnan');
Sj  = 1*median(MX(:,:,6),1,'omitnan');

figure(21);clf;

subplot(121);
plot(12*log2(F1MEAN/440), log10( Vgrv(:)./(co * Spall_Lateral(:))),'o');
hold on; box on; grid on;
plot(12*log2(F1MEAN/440), log10(Vgrv./(co*Sin)) ,'o');
legend('V_{grv}/c_O S_{pall}','V_{grv} / c_o S_{in}');
ylim([-3.5 0]); xlim([-20 23]);

subplot(122);
plot( 12*log2(F1MEAN/440), log10( Vf./(co*Sin)) , 'o');
hold on; box on; grid on;
plot(12*log2(F1MEAN/440), log10(Vf./(co*Sj)), 'o');
legend('V_f / c_o S_{in}','V_f / c_o S_J');
ylim([-2.3 0]); xlim([-20 23]);


% [Vgr/(c * Spall)]

% [Vgr / ( c * Sin)]

% [Vf / (c * Sin)]

% [Vf / ( c * Sj)]


% Characteristic lengths ============

% [L_pall / c]
% 8mm / 340m/s   = 2.3529e-5 s (0.0235 ms)
% [L_in / c]
% ~1mm+ / 340m/s = 0.0029+ ms
% [L_j / c]
% ~1mm/340m/s    = 0.0029 ms

% Geometrical sections comparison =============0
figure();
plot(12*log2(F1MEAN/440));





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Groove and Foot analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% t10 delays: mech, hydro, acoust [[[to keep, 21/01/2025]]]

figure();hold on;
errorbar( median(12*log2(MX(:,:,13)/440),1,'omitnan')', ...
    F1MEAN.*median(abs(MX(:,:,57)-MX(:,:,55) ),1,'omitnan')',...
    std( MX(:,:,13).* abs(   MX(:,:,57)-MX(:,:,55) ) ,'omitnan'),...
    ':v','color',[1,1,1]*0.7);

plot(  median(12*log2(MX(:,:,13)/440),1,'omitnan') , ...
    median(1e3*MX(:,:,51),1,'omitnan'),...
    '-o', 'Color','b');

plot( median(12*log2(MX(:,:,13)/440),1,'omitnan'),...
        median(    1e3*(MX(:,:,55)-MX(:,:,51) ),1,'omitnan'),...
        '-*k');
    
legend('Acoust.Delay/$T_1$','Mech.Delay [ms]','Hydrod.Delay [ms]', 'interpreter','latex','location','northwest');
% ylabel('Delay (lin)','interpreter','latex');
xlabel('$12log_2(f_1/440)$','interpreter','latex'); 
box on; grid on;% ylim([0 2]);
ylim([-3 20]); xlim([-20 23]);

%% PRT10ft and PRT10grv, by tessitura :: with errors [[[Keep 2025/01/30]]]
figure();
errorbar(12*log2(F1MEAN/440), ...
    1e3*median(MX(:,:,52),1,'omitnan'), ...
   std(1e3*MX(:,:,52),1,'omitnan')' ,...
   'marker','d', 'linestyle','none', 'Color', 'k','markerfacecolor','k');

hold on;
errorbar(12*log2(F1MEAN/440), ...
    1e3*median(MX(:,:,56),1,'omitnan'), ...
    std(1e3*MX(:,:,56),1,'omitnan')',...
    'marker','s','linestyle', 'none', 'color', 'k');
    
legend('$PRT^{90}_{grv}$','$PRT^{90}_{ft}$','interpreter','latex','location','NorthEast');
xlabel('$12\times log_2(f_1/440 Hz)$', 'interpreter','latex');
ylabel('[ms]','interpreter','latex');
grid on; box on;
ylim([0 8]);


%% t20grv - t20foot, by tessitura;  06/12/2024[[[[[[06/12/2024 [to keep]]]]]]]

figure();
scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan'), median(1e3*abs(MX(:,:,33)-MX(:,:,34) ),1,'omitnan') , 'dk', 'filled');
ylabel('$t^{20}_{grv}$ - t$^{20}_{foot}$ (lin)[ms]','interpreter','latex');
xlabel('$12log_2(f_1/440)$','interpreter','latex'); 
box on; grid on; ylim([0 2]);

%% PRT20ft and PRT20grv, by tessitura;
figure();
scatter( 12*log2(F1MEAN/440), 1e3*median(MX(:,:,28),1,'omitnan') , 'xk');
hold on;
scatter( 12*log2(F1MEAN/440), 1e3*median(MX(:,:,29),1,'omitnan') , 'ok');
legend('PRT^{20-80}_{grv}','PRT^{20-80}_f');ylim([0 4]);
box on; grid on;


%% t80grv - t80 foot, by tessitura; [[[[[[06/12/2024 [to keep]]]]]]]

figure();
scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan'), ...
    1e3*median(   MX(:,:,29)+MX(:,:,34)-MX(:,:,28)-MX(:,:,33) ,1,'omitnan') ,...
    'dk', 'filled');

ylabel('$t^{80}_{ft}$ - t$^{80}_{grv}$  [ms]','interpreter','latex');
xlabel('$12log_2(f_1/440)$','interpreter','latex'); box on; grid on;
ylim([0 2.4]);
grid on; box on;

%% t80grv - t80foot, by Volume Ratios;  [[[[[[06/12/2024 [to keep]]]]]]]

figure();
scatter( median( MX(:,:,2)./MX(:,:,11),1,'omitnan'), ...
    median(  1e3* abs(MX(:,:,29)+MX(:,:,34)-MX(:,:,28)-MX(:,:,33)) ,1,'omitnan') ,...
    'dk', 'filled');

xlabel('$V_f/V_{grv}$','Interpreter','latex');
ylabel('$t^{80}_f-t^{80}_{grv}$','interpreter','latex');
grid on; box on;
ylim([0 2.5]);

%% t10grv - t10foot, by tessitura;
figure();
plot( 12*log2(F1MEAN/440), 1e3*median(MX(:,:,55)-MX(:,:,51),1,'omitnan') , 'dk','markerfacecolor','k');
grid on; box on;
xlabel('$12\times log_2 (f_1/440 Hz)$','interpreter','latex');
ylabel('$t^{10}_f - t^{10}_{grv}$ [ms]','interpreter','latex');
ylim([0 2.5]);


%% t90ft-t90grv by tessitura
figure();
scatter(  12*log2(F1MEAN/440), 1e3*median(MX(:,:,55)+MX(:,:,56) - MX(:,:,51)-MX(:,:,52),1,'omitnan'), 'kd','filled');
grid on; box on;
xlabel('$12\times log_2(f_1/440 Hz)$','interpreter','latex');
ylabel('$t^{10-90}_{ft}-t^{10-90}_{grv}$ [ms]','interpreter','latex');


%% t90ft-t90grv by Volume ratios
figure();
scatter(  median(MX(:,:,2)./MX(:,:,11),1,'omitnan'), 1e3*median(MX(:,:,55)+MX(:,:,56) - MX(:,:,51)-MX(:,:,52), 1, 'omitnan'), 'kd','filled' );
grid on; box on;

xlabel('$V_f/V_{grv}$','interpreter','latex');
ylabel('$t^{90}_{ft} - t^{90}_{grv}$ [ms]','interpreter','latex');





%% t20_foot normalized by tau_f = Vf/(c_o * S_in) [HYDRODYNAMIC DELAY]

figure();
scatter( 12*log2(F1MEAN/440),median( co * MX(:,:,34).*MX(:,:,5)./MX(:,:,2) ,1,'omitnan'), 'b', 'filled');
ylabel('$ t^{20}_f / \tau_f $ with $\tau_f = V_f / c_o  \mathcal{S}_{in}$','interpreter','latex');
xlabel('$12\times log_2(f_1/440)$','interpreter','latex'); box on; grid on;
title('Hydrodynamic delay','interpreter','latex');
ylim([0 1.3]);

%% t20_foot normalized by tau_f = Vf/(u_in * S_in) [HYDRODYNAMIC DELAY]

tau_f = MX(:,:,2)./(MX(:,:,5).*sqrt( 2/rho*(MX(:,:,20)-MX(:,:,21))) );

figure();
scatter( 12*log2(F1MEAN/440), median(MX(:,:,34)./tau_f,1,'omitnan') , 'b', 'filled');
ylabel('$ t^{20}_f / \tau_f $ with $\tau_f = V_f / u_{in}  \mathcal{S}_{in}$','interpreter','latex');
xlabel('$12\times log_2(f_1/440)$','interpreter','latex'); box on; grid on;
title('Hydrodynamic delay','interpreter','latex');
ylim([0 0.2]);

%% t20_foot normalized by tau_f = Vf/(u_j * S_j) [HYDRODYNAMIC DELAY]

tau_f = MX(:,:,2)./(MX(:,:,6).*sqrt( 2/rho*(MX(:,:,21)-0*MX(:,:,21))) );

figure();
scatter( 12*log2(F1MEAN/440), median(MX(:,:,34)./tau_f,1,'omitnan') , 'b', 'filled');
xlabel('Tessitura [semitones]');
ylabel('$ t^{20}_f / \tau_f $ with $\tau_f = V_f / u_{j}  \mathcal{S}_{j}$','interpreter','latex');
grid on; box on;
ylim([0 0.055]);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acoustic analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%% PRTrad/T1 (w.r.t tessitura) [OK]
figure();
errorbar(    median(12*log2(MX(:,:,13)/440),1,'omitnan'), ...
    median(MX(:,:,30).*MX(:,:,13), 1, 'omitnan'),...
    std(MX(:,:,30).*MX(:,:,13), 'omitnan'), ...
    'linestyle','none', 'marker', 'o', 'color', 'k', 'linewidth', 1); 

xlabel('$12 \times log_2(f_1/f_{ref})$','interpreter','latex');
ylabel('$ PRT_{rad} /T_1$','interpreter','latex');
box on; grid on;



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


%% PRTfoot vs. PRTgroove [OK]

figure();

scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan'), median(MX(:,:,29)./MX(:,:,28),1,'omitnan'), 'k', 'filled');
ylabel('PRT_{foot}/PRT_{groove} [n.u.]');
grid on; box on;
ylim([0 1.3]);
yyaxis right;
scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan'), median(1e3*(MX(:,:,29)-MX(:,:,28)),1,'omitnan'), ...
    'markeredgecolor', "#D95319", 'MarkerFaceColor',"#D95319");
ylabel('PRT_{foot} - PRT_{groove} [ms]');
ylim([-1.5 2]);
xlabel('tessitura [semitones]');


%% t20 differences by tessitura 

%>>> tau = Vf/(c_o * S_in)
tau = median(MX(:,:,2)./(co.*MX(:,:,5)), 1, 'omitnan');
tau = 1;
figure();
errorbar(12*log2(F1MEAN/440),...
     1e3*median(MX(:,:,34)-MX(:,:,33),1,'omitnan')./tau,...
    1e3*std(  (MX(:,:,34)-MX(:,:,33)),1,'omitnan' ) ./ tau,...
    'linestyle','none','color','k','marker', 'o', 'linewidth',1);
ylabel('$t^{20}_{f}-t^{20}_{grv}$','interpreter','latex');
xlabel('Tessitura [semitones]','interpreter','latex');
title('t20 differences');
box on; grid on;
 ylim([0 2]);

%% PRT ratios by tessitura 
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
    1e3*median(MX(:,:,29),1,'omitnan'),...
    1e3*std(MX(:,:,29)),...
    'linestyle','none','color','k','marker', 'o','linewidth',1);
ylim([0 4]);
xlabel('Vf'); ylabel('PRT_f [ms]');




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transient spectrum study
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% max(a2/a1) during transient [not bad]

% GVC: I have used butter(4), filtfilt() of envel_first and envel_second
% and then search max within t20_f and t80_f+50*PRT.

mask = [1:22]; %mask(7)=[]; % >>> If you want to remove a single dataset


figure();
errorbar(12*log2(F1MEAN/440), median(MX(:,:,48),1,'omitnan'), std(MX(:,:,48),1,'omitnan'));
grid on;
% ylabel('NOT THIS: trapez(dt, a_2(t) - a_1(t) )/T1  ,  t=[t^{20}_{p}, t^{80}_{p}] [n.u.]'); xlabel('tessitura semitones');
xlabel('tessitura');
grid on; box on;
title('Integrated area comprised between a_2(t) and a_1(t) in the interval (t^{20}_p, t^{80}_p) norm. by oscill. period');


%% Transient spectra differences integrated in area between a_3(t) and a_1 for t=[t20_p, t80_fp [not bad]

figure();
errorbar( 12*log2(MX(:,:,13)/440),    median(MX(:,:,39).*MX(:,:,13),1,'omitnan'), std(MX(:,:,39).*MX(:,:,13),1,'omitnan'),'kd' );
% hold on;
% plot( 12*log2(F1MEAN/440), median( MX(:,:,39).*MX(:,:,13), 1, 'omitnan')  , '-ro');

ylabel('trapez(dt, a_3(t) - a_1(t) )/T_1, t=[t^{20}_p, t^{80}_p)] [n.u.]'); xlabel('tessitura [semitones]');

grid on; box on;
title('Integrated area comprised between a_3(t) and a_1(t) in the interval (t^{20}_p, t^{80}_p) norm. by oscill. period');
ylim([-5 20]);


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

%%
clear MXresh;
MXresh = [];
for idx = 1 : size(MX,3)
% for idx = 14 : size(MX,3)    
    tmp =  MX(:,:,idx);
    tmp = tmp(:);
    tmp(isnan(tmp)) = [];
    tmp(isinf(tmp)) = [];
    tmp = squeeze(tmp);
    tmp = tmp(1:908);
   MXresh =  [MXresh, tmp];  
end
MXresh = (MXresh-mean(MXresh,1))./std(MXresh,1);

% linear pca()
[coeff,score,latent,tsquared,explained,mu] = pca(MXresh);
% log pca()
MXresh = MXresh - min(MXresh) + 1e-6;
MXresh_log = log10(MXresh);
MXresh_log = (MXresh_log-mean(MXresh_log,1))./std(MXresh_log);
[coeff_log,score_log,latent,tsquared,explained_log,mu] = pca(MXresh_log);

figure(72);clf;
subplot(1,4,[1:3]);
heatmap( (coeff), 'Colormap', jet(4), 'ColorbarVisible', true);
title('PCA lin'); xlabel('P components');ylabel('Variables');
colormap(bluewhitered);clim([-1,1]);
subplot(1,4,4);
pareto(explained); grid on; 

figure(73);clf;
subplot(1,4,[1:3]);
heatmap( (coeff_log), 'Colormap', jet(4), 'ColorbarVisible', true);
title('PCA log'); xlabel('P components');ylabel('Variables');
colormap(bluewhitered);clim([-1,1]);
subplot(1,4,4);
pareto(explained_log); grid on;


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

maskpca = 1:50;
DataToWorkWith = MXresh(:,maskpca);


% tmp = (tmp-mean(tmp,1,'omitnan'))./var(tmp,'omitnan');
tmp = (tmp-mean(tmp,1,'omitnan'))./std(tmp,'omitnan');

% tmp = exp(tmp);
%%

tmp = MXresh;

tmpc = corrcoef(tmp);
tmpc(abs(tmpc)<0.0001) = 0;



namevarsall = {'Lp','Vf','PWdth','TNHD','INLET','SJET','h','H','Wm','Dp',...
    'Vgrv','Qfact','f1','theta','Qpall2grv','Qgr2ft','Qjet','Rmplssg',...
    'P{o} pall','P{o}grv','P{o}foot','P{o}rad',...
    'P{o} grv-pall','P{o}foot-grv','P{o}rad-foot',...
    'Beta','nu',...
    'PRTgrv','PRTfoot','PRTrad',...
    'PRTfoot/grv','PRTrad/foot',...
    't{20}grv','t{20}foot','t{20}rad',...
    't20(foot-grv)','t20(rad-foot)',...
    'Area1','Area2','Sin/Spal','Sjet/Sin','MaxKeyVel',...
    'a2maxOverA1simult','a2maxOvera1Targ','a2maxOvera2Targ','DeltaPfoot2mouthAta2max','a2maxVec','a2maxOverA1smooth',...
    'gofr2','LateralChkSec'};





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

% Jordi Tur 2024/12/13:
% y ~ p1 + p2 + ... pk
% log(y) ~ p1 + p2 + ... pk
% log(y) = b1*p1 + b2*p2 + ... bk*pk
% y = exp( b1*p1 + b2*p2 + ... bk*pk ) =
%   = exp( b1*p1 ) exp( b2*p2 ) * ... * exp( bk*pk )
% 
%   bi <- exp(bi)
%   exp(bi*(xi+0.1)) / exp(bi*xi)


maskpca = [14:56];

DataToWorkWith = MXresh(:,maskpca);


olo = (DataToWorkWith-mean(DataToWorkWith,1,'omitnan'))./std(DataToWorkWith,'omitnan');
olo = olo+1.05;
olo = log10(olo);

%%

%[coeff,score,latent,tsquared,explained] = pca(DataToWorkWith, 'VariableWeights',wei);
[coeff,score,~,~,explained] = pca(MX(1:40,:,:));
figure(3);clf; plot(explained,'o-');
LBL = append(namevarsall([13:50]));


figure(4);clf;
%biplot(coeff(:,1:3),'scores',score(:,1:3),'VarLabels',LBL);
biplot(coeff(:,1:3),'scores',score(:,1:3));


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
