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


rho = 1.2; % Air density [kg/m^3]
co  = 340; % Acoustic propagation speed [m/s]
co2 = co^2; 
P0  = 820; % Pallet box pressure [Pa = kg m^-1 s^-2]

Vgrv = median(MX(:,:,11),1,'omitnan');
Vf   = median(MX(:,:,2),1, 'omitnan');

Spall_Slot    = median(MX(:,:,3),1,'omitnan')*0.1298; % Perforated rectangles on the plate, with the same width as the groove channel
Spall_Lateral = PalletValveStrokeArea(maskpipes); % At maximum aperture of valve, adding areas of a rectangle and two triangles

fax = 12*log2(F1MEAN/440);


%%

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %       Geometry
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All geometry vars w.r.t. f1/440
FSZ = 17;
figure(); hold on;

for idx = 1 :12
    tmp_data = log10( median(MX(:,:,idx)) );
    val_absc = interp1(freqlogax, tmp_data, 0);
    plot(freqlogax, tmp_data-val_absc ,'-') ;    
end
grid on; box on;
legend(namevarsall);

if 0%figure(); 
% for idx=1:12
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

%% S_palletvalve, S_slot, S_tonehole, S_in, S_j 
figure(); hold on; grid on; box on;

scatter(freqlogax, median(1e6*MX(:,:,50),1,'omitnan'), 'b', 'filled');       % mean = 8.4 cm^2 % pallet valve stroke
scatter(freqlogax, median(1e6*MX(:,:,3)*0.1298,1,'omitnan'), 'r', 'filled'); % mean = 16.3 cm^2 % pallet slot window
scatter(freqlogax, median(1e6*pi*(0.5*MX(:,:,4)).^2,1,'omitnan'),'k','filled' ); % tone hole area
scatter(freqlogax, median(1e6*MX(:,:,5),1,'omitnan'), 'g', 'filled'); % foot inlet
scatter(freqlogax, median(1e6*MX(:,:,6),1,'omitnan'), 'm', 'filled'); % foot outlet
ylabel('mm^2'); xlabel('tessitura [semitones]');
ax=gca; ax.YScale = 'log';
legend('Lateral PalletValve','Pallet Slot','Tone Hole','Sin','Sj');

% Sections
            % % Geometrical sections comparison =============0
            % figure(22);clf;
            % plot(12*log2(F1MEAN/440), log10(Spall_Lateral)   ,'-*');
            % hold on;
            % plot(12*log2(F1MEAN/440), log10(Sin), '-*');
            % plot(12*log2(F1MEAN/440), log10(Sj), '-*');
            % legend('S_{pall,lat} (geom)','S_{in}    (geom)','S_j    (geom)');
            % box on; grid on; 
            % ylim([-6 0]); ylabel('log10');


% Sj / Sin

figure();
scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan'), median((MX(:,:,6)./MX(:,:,5)),1,'omitnan'), 'kd', 'filled');
ylabel('$\mathcal{S}_j / \mathcal{S}_{in} \ (lin)$ [n.u.]', 'interpreter','latex');
xlabel('$12log_2(f_1/440)$','interpreter','latex'); box on; grid on; ylim([0 1.5]);



%%                          % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %       Steady-State analysis
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target pressures study
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ptarg differences all wrt to reservoir, in SS grv and rsv have the same pressure [OK,final][2025/01/30]
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

% Theta [OK]
figure();

scatter( median(12*log2(MX(:,:,13)/440),1,'omitnan'), median(MX(:,:,14),1,'omitnan'), 'dk', 'filled');
grid on; box on;
xlabel('$12 \times log_2 (f_1 /  440 Hz) $', 'interpreter','latex');
ylabel('$\theta = u_j/f_1 W_m$','interpreter','latex'); ylim([0 12]);

% Flow conservation for pallet-to-groove

Spall_eff = median(MX(:,:,6),1,'omitnan').*sqrt(median(MX(:,:,21),1,'omitnan'))./sqrt( median(MX(:,:,19),1,'omitnan')-median(MX(:,:,20),1,'omitnan'));

figure();
Pallet_pfit = polyfit(fax,log10(Spall_eff),1);
plot(fax, log10(Spall_eff), '-o');
ylabel('Spall eff (log10)');
grid on;
hold on;

querypoints = [-20:0.1:23];
plot(querypoints, polyval(Pallet_pfit, querypoints), '--k');
text(0.1,-4.4,sprintf('$y=%1.5f x + %1.4f$',Pallet_pfit(1),Pallet_pfit(2) ), 'interpreter','latex');
xlabel('12log2(f1/440)');

figure();
Sratio = median(MX(:,:,6),1,'omitnan')./Spall_Lateral';

plot(Sratio, Spall_eff, '-o');
xlabel('Sj/Spall');
ylabel('Spall Effective');
grid on;

Pallet_pfit2 = polyfit(Sratio,Spall_eff,1);
querypoints = [0:1e-4:0.08];
hold on;
plot(querypoints, polyval(Pallet_pfit2, querypoints), '--k');
text(0.02,0.7e-4,sprintf('$y=%1.5f x + %1.4f$',Pallet_pfit(1),Pallet_pfit(2) ), 'interpreter','latex');


% Foot flow conservation and Gamma function [OK][Keep, 2025/01/30, plot1/2 venacontracta]
      
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
% Pfoot target versus Sj/Sin [ok][keep][2025/01/30][plot2/2, venacontracta]

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

 % Additional hypothesis on why we don't see the flow consrevation:
 % the foot is a big volume with a small channel at the entrance and output
 % so we have a jet at the entrance dissipated by turbulence
 % and a second at the output....? so?

%% P_0 by tessitura
figure;
errorbar(fax,...
    mean(MX(:,:,19),1,'omitnan'),...
    std(MX(:,:,19),1,'omitnan'),...
    '-v');
xlabel('12log_2(f_1/440)');ylabel('P_0 [Pa]'); grid on; box on;ylim([0 900]);

%%                      % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %       Transient analysis
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===     Characteristic times: repetition of below's (including A, B, C, D)
fax = 12*log2(F1MEAN/440);
LW = 1.0;

l_pall = 5e-3;
l_in   = 1e-3;
l_j    = 1e-3;

Sgrv_eff = median(MX(:,:,3),1,'omitnan')'*0.05;


Vgrv = median(MX(:,:,11),1,'omitnan'); Vgrv = Vgrv(:);
Vf   = median(MX(:,:,2),1, 'omitnan'); Vf = Vf(:);

Spall_max = PalletValveStrokeArea(maskpipes);
    Spall_max = Spall_max(:); 
Sin           = MX(:,:,5);
    Sin = median(Sin, 1, 'omitnan'); Sin = Sin(:);
Sj            = 1*MX(:,:,6);
    Sj = median(Sj, 1, 'omitnan'); Sj = Sj(:);


        % Spall_max_eff = 1.*Spall_max;
Spall_max_eff = 10.^(  -0.03732*12*log2(F1MEAN/440) - 4.6202 ); % fitted to computed 
Sin_eff = Sj(:).*sqrt(median(MX(:,:,21),1,'omitnan')'./(median(MX(:,:,19),1,'omitnan')- median(MX(:,:,21),1,'omitnan') )'         );
Sj_eff = 1.*Sj;

tau_pall_L = l_pall*sqrt(rho/P0)*ones(length(Sin),1);
tau_in_L   = l_in  *sqrt(rho/P0)*ones(length(Sin),1);
tau_j_L    = l_j   *sqrt(rho/P0)*ones(length(Sin),1);

tau_pall_V = Vgrv./(Spall_max_eff*co2)*sqrt(P0/rho);
tau_in_V   = Vgrv./(Sin_eff*co2)*sqrt(P0/rho);
tau_j_V    = Vgrv./(Sj_eff*co2)*sqrt(P0/rho);

Vf_over_Vgrv = Vf./Vgrv;

one_over_Amax = Spall_max_eff*co2./Vgrv*sqrt(rho/P0);
one_over_B    = Sin_eff*co2./Vgrv*sqrt(rho/P0);
% C = Sin_eff*co2./Vgrv*sqrt(rho/P0);
% D = Sj_eff*co2./Vgrv*sqrt(rho/P0);
one_over_C = Sin_eff*co2./Vf*sqrt(rho/P0);
one_over_D = Sj_eff*co2./Vf*sqrt(rho/P0);
sigMa = Spall_max_eff(:)./Sgrv_eff(:);


figure(24); clf;
plot(fax, median(MX(:,:,52),1,'omitnan') ,'-o', 'linewidth',LW);
hold on;
plot(fax, median(MX(:,:,56),1,'omitnan') ,'-o', 'linewidth',LW);
% plot(fax, median(MX(:,:,58),1,'omitnan') , '-o', 'linewidth',LW);
%
plot(fax, tau_pall_L, '--s', 'linewidth',LW);
plot(fax, tau_in_L , '--s', 'linewidth',LW);
plot(fax, tau_j_L,'--s', 'linewidth',LW);
%----------
A = 1./one_over_Amax;
B = 1./one_over_B;
C = 1./one_over_C;
D = 1./one_over_D;
plot(fax, A,'-v', 'linewidth',LW);
plot(fax, B,'-v', 'linewidth',LW);
plot(fax, C,'-v', 'linewidth',LW);
plot(fax, D,'-v', 'linewidth',LW);
% plot(fax, one_over_Amax,'-v', 'linewidth',LW);
% plot(fax, one_over_B,'-v', 'linewidth',LW);
% % plot(fax, Vf./(Vgrv.*C),'-v', 'linewidth',LW;
% % plot(fax, Vf./(Vgrv.*D),'-v', 'linewidth',LW);
% plot(fax, one_over_C,'-v', 'linewidth',LW);
% plot(fax, one_over_D,'-v', 'linewidth',LW);
%
legend('PRT_{grv}','PRT_f','tau grv (L_{grv})','tau in (L_{in})','tau jet (L_{j})','A_{max}','B','C','D', 'fontsize',12);
ax=gca; ax.YScale = 'log'; grid on;
xlabel('12log_2(f_1/440)');


%%
% Fig4 A
figure();
plot(Sin./PalletValveStrokeArea(maskpipes),...
    A./B,...
    'dk', 'markerfacecolor','k');
ylabel('A/B');
xlabel('S^{geo}_{in}/S^{geo}_{pall}');
grid on; ylim([0 0.55]);xlim([0 0.07]);


%%
% Fig4B
figure();
plot(Sj./Sin,...
    C./D,...
    'sk','markerfacecolor','k');
grid on; 
ylim([0 1.8]);ylabel('C/D');
xlim([0 1.4]);xlabel('S^{geo}_j/S^{geo}_{in}');

                   
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulated vs Measured ABCDsigma params 

A = 1./one_over_Amax;
B = 1./one_over_B;
C = 1./one_over_C;
D = 1./one_over_D;
% sigMa = sigMa;

A = A(:); B = B(:); C = C(:); D = D(:); sigMa = sigMa(:);
PRTgrv = median(MX(:,:,52),1,'omitnan');
PRTgrv = PRTgrv(:);
PRTf = median(MX(:,:,56), 1, 'omitnan');
PRTf = PRTf(:);

pf_over_pgrv = median(MX(:,:,21),1,'omitnan')./median(MX(:,:,20),1,'omitnan');
pf_over_pgrv = pf_over_pgrv(:);
%
figure();
subplot(121);
plot(A, 1e3*PRTgrv,'d'); hold on;  plot(A, 1e3*PRTf,'s'); 
xlabel('A [s]'); ylabel('[ms]'); legend('PRT_{grv}','PRT_{f}'); title('PRT'); grid on; box on;
subplot(122);
plot(A, pf_over_pgrv,'v'); xlabel('A [s]');title('Pf/P_{grv}'); grid on; box on; ylim([0 1]);
%
figure();
subplot(121);
plot(B, 1e3*PRTgrv,'d'); hold on;  plot(B, 1e3*PRTf,'s'); 
xlabel('B [s]'); ylabel('[ms]'); legend('PRT_{grv}','PRT_{f}'); title('PRT'); grid on; box on;
subplot(122);
plot(B, pf_over_pgrv,'v'); xlabel('B [s]');title('Pf/P_{grv}'); grid on; box on; ylim([0 1]);
%
figure();
subplot(121);
plot(C, 1e3*PRTgrv,'d'); hold on;  plot(C, 1e3*PRTf,'s'); 
xlabel('C [s]'); ylabel('[ms]'); legend('PRT_{grv}','PRT_{f}'); title('PRT'); grid on; box on;ylim([0 8]);
subplot(122);
plot(log10(C), log10(pf_over_pgrv),'v'); xlabel('C [s]');title('Pf/P_{grv}'); grid on; box on; %ylim([0 1]);
%
figure();
subplot(121);
plot(D, 1e3*PRTgrv,'d'); hold on;  plot(D, 1e3*PRTf,'s'); 
xlabel('D [s]'); ylabel('[ms]'); legend('PRT_{grv}','PRT_{f}'); title('PRT'); grid on; box on;ylim([0 8]);
subplot(122);
plot(log10(D), log10(pf_over_pgrv),'v'); xlabel('D [s]');title('Pf/P_{grv}'); grid on; box on;% ylim([0 1]);
%
figure();
subplot(121);
plot(sigMa, 1e3*PRTgrv,'d'); hold on;  plot(sigMa, 1e3*PRTf,'s'); 
xlabel('\Sigma'); ylabel('[ms]'); legend('PRT_{grv}','PRT_{f}'); title('PRT'); grid on; box on;
subplot(122);
plot(sigMa, pf_over_pgrv,'v'); xlabel('\Sigma');title('Pf/P_{grv}'); grid on; box on; ylim([0 1]);





%% ===     Characteristic times: length, volume, and volume ratio vs. PRT's [potentially yes]
fax = 12*log2(F1MEAN/440);
LW = 1.0;

l_pall = 5e-3;
l_in   = 1e-3;
l_j    = 1e-3;

Vgrv = median(MX(:,:,11),1,'omitnan'); Vgrv = Vgrv(:);
Vf   = median(MX(:,:,2),1, 'omitnan'); Vf = Vf(:);

Spall_max = PalletValveStrokeArea(maskpipes);
    Spall_max = Spall_max(:); Spall_max = Spall_max(:);
Sin           = MX(:,:,5);
    Sin = median(Sin, 1, 'omitnan'); Sin = Sin(:);
Sj            = 1*MX(:,:,6);
    Sj = median(Sj, 1, 'omitnan'); Sj = Sj(:);


% Spall_max_eff = 1.*Spall_max;
Spall_max_eff = 10.^(  -0.03732*12*log2(F1MEAN/440) - 4.6202 ); % fitted to computed 
Sin_eff = Sj(:).*sqrt(median(MX(:,:,21),1,'omitnan')'./(median(MX(:,:,19),1,'omitnan')- median(MX(:,:,21),1,'omitnan') )'         );
Sj_eff = 1.*Sj;

tau_pall_L = l_pall*sqrt(rho/P0)*ones(length(Sin),1);
tau_in_L   = l_in  *sqrt(rho/P0)*ones(length(Sin),1);
tau_j_L    = l_j   *sqrt(rho/P0)*ones(length(Sin),1);

tau_pall_V = Vgrv./(Spall_max_eff*co2)*sqrt(P0/rho);
tau_in_V   = Vgrv./(Sin_eff*co2)*sqrt(P0/rho);
tau_j_V    = Vgrv./(Sj_eff*co2)*sqrt(P0/rho);

Vf_over_Vgrv = Vf./Vgrv;


figure(26); clf;
plot(fax, median(MX(:,:,52),1,'omitnan') ,'-o', 'linewidth',LW);
hold on;
plot(fax, median(MX(:,:,56),1,'omitnan') ,'-o', 'linewidth',LW);
plot(fax, median(MX(:,:,58),1,'omitnan') , '-o', 'linewidth',LW);
%
plot(fax, tau_pall_L, '--s', 'linewidth',LW);
plot(fax, tau_in_L , '--s', 'linewidth',LW);
plot(fax, tau_j_L,'--s', 'linewidth',LW);
%
plot(fax, 1.*tau_pall_V,'-.*', 'linewidth',LW);
plot(fax, 1.*tau_in_V, '-.*', 'linewidth',LW);
plot(fax, 1.*tau_j_V, '-.*', 'linewidth',LW);
%
plot(fax, 1.*Vf_over_Vgrv,'-vk', 'linewidth',LW);
%
legend('PRT_{grv}','PRT_f','PRT_m','tau grv (L)','tau in (L)','tau jet (L)','tau grv (V)','tau f (V)','tau jet (V)','Vf/Vgrv', 'fontsize',12);
ax=gca; ax.YScale = 'log'; grid on;
xlabel('12log_2(f_1/440)');

% =========================
figure(27); clf;plot(fax, median(MX(:,:,52),1,'omitnan') ,'-o', 'linewidth',LW);
hold on;
plot(fax, median(MX(:,:,56),1,'omitnan') ,'-o', 'linewidth',LW);
plot(fax, median(MX(:,:,58),1,'omitnan') , '-o', 'linewidth',LW);
%
plot(fax, tau_pall_L, '--s', 'linewidth',LW);
plot(fax, tau_in_L , '--s', 'linewidth',LW);
plot(fax, tau_j_L,'--s', 'linewidth',LW);
%
plot(fax, 1./tau_pall_V,'-.*', 'linewidth',LW);
plot(fax, 1./tau_in_V, '-.*', 'linewidth',LW);
plot(fax, 1./tau_j_V, '-.*', 'linewidth',LW);
%
plot(fax, 1.*Vf_over_Vgrv,'-vk', 'linewidth',LW);
%
legend('PRT_{grv}','PRT_f','PRT_m','tau grv (L)','tau in (L)','tau jet (L)','tau grv (V)','tau f (V)','tau jet (V)','Vf/Vgrv', 'fontsize',12);
ax=gca; ax.YScale = 'log'; grid on;
xlabel('12log_2(f_1/440)');



%% testing A, B, C, D eigs

M = zeros(2,2,22);
eigmax = zeros(22,1);
eigmin = zeros(22,1);
for idx = 1:22
    M(:,:,idx) = [A(idx), B(idx); C(idx), D(idx)];
    tmp = eigs(M(:,:,idx));
    eigmax(idx) = tmp(1);
    eigmin(idx) = tmp(2);
end

figure();
% plot(A*1e3, 1e3*eigmax);
plot(1e3*eigmax, A*1e3 , 'o');
hold on;
plot(PRTf*1e3, A*1e3, 'o');
plot(PRTf*1e3, eigmax*1e3, 'o');

plot([0 20],[0 20],'-k');
grid on; axis equal;
%%
figure();
subplot(211);
plot(A*1e3);
hold on;
plot(eigmax*1e3);
grid on;
subplot(212);
plot(A./eigmax); grid on;
%%
figure(31); hold on;
plot(eigmax);
plot(A);
plot(B);
plot(C);
plot(D);
plot(PRTf);
ax=gca; ax.YScale = 'log'; grid on; box on;
legend('eigmax','A','B','C','D', 'PRTf');


            %% PRTgrv vs Vgrv, PRTf vs Vf [NO]
            figure();
            %
            subplot(121);
            plot(median(MX(:,:,11),1,'omitnan'),median(MX(:,:,52),1,'omitnan') ,'ko','markerfacecolor','k');
            xlabel('Vgrv'); ylabel('PRTgrv');
            ylim([0 7.5e-3]);
            xlim([0 4.5e-4]);
            grid on;
            %
            subplot(122);
            plot(median(MX(:,:,2),1,'omitnan'), median(MX(:,:,56),1,'omitnan') ,'ko','markerfacecolor','k');
            xlabel('Vf'); ylabel('PRT f');
            ylim([0 7e-3]);
            xlim([0 3.5e-4]);
            grid on;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Groove and Foot analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% t90 groove, t90 foot [just for inspection]
figure();

errorbar(12*log2(F1MEAN/440),...
    median( 1e3*(MX(:,:,51)+MX(:,:,52)),1,'omitnan' ),...
    std( 1e3*(MX(:,:,51)+MX(:,:,52)) ,1,'omitnan') ) ;

hold on;

errorbar(12*log2(F1MEAN/440),...
    median( 1e3*(MX(:,:,55)+MX(:,:,56)), 1, 'omitnan') ,...
    std(  1e3*(MX(:,:,55)+MX(:,:,56)) , 1, 'omitnan') );

grid on; box on; 
ylim([0 16]);
legend('t^{90}_{grv}','t^{90}_{ft}');

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

%% t90ft-t90grv by tessitura
figure();
scatter(  12*log2(F1MEAN/440), 1e3*median(MX(:,:,55)+MX(:,:,56) - MX(:,:,51)-MX(:,:,52),1,'omitnan'), 'kd','filled');
grid on; box on;
xlabel('$12\times log_2(f_1/440 Hz)$','interpreter','latex');
ylabel('$t^{10-90}_{ft}-t^{10-90}_{grv}$ [ms]','interpreter','latex');


%% t90ft-t90grv by Volume ratios
figure();
scatter(  1./median(MX(:,:,2)./MX(:,:,11),1,'omitnan'), 1e3*median(MX(:,:,55)+MX(:,:,56) - MX(:,:,51)-MX(:,:,52), 1, 'omitnan'), 'kd','filled' );
grid on; box on;

xlabel('$V_f/V_{grv}$','interpreter','latex');
ylabel('$t^{90}_{ft} - t^{90}_{grv}$ [ms]','interpreter','latex');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acoustic analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% empty

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



                %% a2max over a2-target [OK]
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
% Study of Beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% beta vs PRTfoot [OK]
figure();
scatter( 1./MX(:,:,26), 1e3*MX(:,:,29)  , 'b','filled');
xlabel('$\beta^{-1} [s]$', 'interpreter','latex');
ylabel('PRT$_{foot}$ [ms]', 'Interpreter','latex');
box on; grid on;
xlim([0 2.2e-3]);
ylim([0 4.5]);

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

%% Beta x T1 vs freq errorbar()


% ...

%% ??????????????????????????
figure();
scatter( 12*log2( MX(:,:,13)/440 ), log10(MX(:,:,17)) ,'filled');
ylim([-5 0]);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Old vena-contracta investigations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

