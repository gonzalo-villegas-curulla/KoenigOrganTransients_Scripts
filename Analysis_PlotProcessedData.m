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


rho = 1.2;
co = 340;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Geometry
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sj / Sin [OK]

figure();
scatter( 12*log2(MX(:,:,13)/440), log10(MX(:,:,6)./MX(:,:,5) ), 'b', 'filled');
ylabel('$\mathcal{S}_j / \mathcal{S}_{in}$ (log10)');
xlabel('$12log_2(f_1/440)$'); box on; grid on;

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Steady-State analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target pressures study
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Ptarg differences [OK]
figure();

scatter(12*log2(MX(:,:,13)/440), abs(MX(:,:,19)-MX(:,:,20)) ,'blue','filled');hold on;
% figure();
scatter(12*log2(MX(:,:,13)/440), abs(MX(:,:,20)-MX(:,:,21)) ,'red','filled');hold on;
% figure();
% scatter(12*log2(MX(:,:,13)/440), abs(MX(:,:,21)-0 ), 'green','filled');

xlabel('$12 \times log_2(f_1/f_{ref})$','interpreter','latex');
ylabel('Target pressure differences','interpreter','latex');
title('Ptarg differences: blue=reserv2grv, red=grv2foot,');
grid on; box on;

%% Ptarg ratios [OK]
figure();

scatter(12*log2(MX(:,:,13)/440), abs(MX(:,:,20)./MX(:,:,19)) ,'blue','filled');hold on;
% figure();
scatter(12*log2(MX(:,:,13)/440), abs(MX(:,:,21)./MX(:,:,20)) ,'red','filled');hold on;
% figure();
% scatter(12*log2(MX(:,:,13)/440), abs(MX(:,:,22)./MX(:,:,21)), 'green','filled');

xlabel('$12 \times log_2(f_1/f_{ref})$','interpreter','latex');
ylabel('Target pressure differences','interpreter','latex');
title('Ptarg ratios: blue=reserv2grv, red=grv2foot,');
grid on; box on;




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


%% t20foot - t20groove [groove delay]

figure();
scatter( 12*log2(MX(:,:,13)/440) , 1e3*( MX(:,:,34) - MX(:,:,33)), 'b', 'filled');
xlabel('tessitura semitones');
ylabel('t20foot - t20groove [ms]');
box on; grid on; ylim([0 2.2]);

yyaxis right;
scatter( 12*log2(MX(:,:,13)/440), 1e3*MX(:,:,33),'color', "#D95319", 'MarkerFaceColor',	"#D95319");
ylabel('t20groove [ms]');
ylim([0 10]);
title('groove-foot delay');


%% PRTfoot vs. PRTgroove [  ]

figure();

scatter( 12*log2(MX(:,:,13)/440), MX(:,:,29)./MX(:,:,28), 'b', 'filled');
ylabel('PRTfoot/PRTgroove');
grid on; box on;
ylim([0 1.3]);
yyaxis right;
scatter( 12*log2(MX(:,:,13)/440), 1e3*(MX(:,:,29)-MX(:,:,28)), 'color', "#D95319", 'MarkerFaceColor',"#D95319");
ylabel('PRTfoot - PRTgroove [ms]');
ylim([-1.5 2]);
xlabel('tessitura [semitones]');

%% Vf vs Vgroove [  ]

figure();

scatter( 12*log2(MX(:,:,13)/440), log10(MX(:,:,11)./MX(:,:,2)), 'b', 'filled');
box on; grid on;
ylabel('Vfoot / Vgroove (log10)');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transient spectrum study
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
hold on;
plot( 12*log2(F1MEAN/440), median( MX(:,:,39).*MX(:,:,13), 1, 'omitnan')  , '-ro');

ylabel('T1 * trapez(dt, a_3 - a_1), (t20_p, t80_p) '); xlabel('tessitura semitones');

grid on; box on;
title('Integrated area comprised between a3 and a1 in the interval (t^{20}_p, t^{80}_p) normalized by oscill. period');
ylim([-5 40]);


%% Fig6, CFA Ernoult2016

mask = [1:22]; mask([1,2])=[];

Uj_at_a2max =  sqrt(MX(:,mask,46)*2/1.2);
theta_a2max = Uj_at_a2max./(MX(:,mask,9).*MX(:,mask,13));
a2_nondim   = MX(:,mask,47)./(rho*co*Uj_at_a2max);

figure();
% scatter( 12*log2(MX(:,:,13)/440),a2_nondim);
scatter(theta_a2max, a2_nondim,'filled');
xlabel('theta at t=a2max');
ylabel('a2max / rho c_o Uj at a2max');
grid on; box on;
xlim([6 11]);ylim([0 2.2e-3]);



%% a2max over a2 target
figure();
scatter( 12*log2(MX(:,:,13)/440), MX(:,:,45), 'b', 'filled');
xlabel('freq in semitones');
ylabel('a2 max over a2 target');
grid on; box on;ylim([0 5]);

%% a2max/a1(@a2max) (mouth rad) [OK]
figure();
scatter(12*log2(MX(:,:,13)/440), (MX(:,:,43)), 'b','filled' );
xlabel({'Sample pipe xticks','12log2(f1/440) spacing'},'interpreter','latex');
ylabel('a2max/a1(@a2max) (log$_{10}$)','interpreter','latex'); 
box on; grid on;

figure();
funh = boxplot( log10(abs(MX(:,:,43))) , 'Positions',12*log2(F1MEAN/440), 'Widths',1);
ylabel('a2max / a1 simult. (log10)');



%% a2max/a1target (mouth rad) [OK]

figure();
scatter(12*log2(MX(:,:,13)/440), (MX(:,:,44)), 'b','filled' );
xlabel('$12*log_2(f_1/f_{ref})$','interpreter','latex');
ylabel('a2max/a1target (lin)','interpreter','latex'); 
box on; grid on; ylim([0 2.2]);

figure();boxplot( (MX(:,:,44)), 'Positions',12*log2(F1MEAN/440), 'Widths',1);
xlabel({'sample pipe xticks','12 log2(f1/440) spacing'},'interpreter','latex');
ylabel('a2max/a1target (lin)','interpreter','latex');
 ylim([0 2.2]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pmouth as per massage to equations 9-10-11 and alpha_vc = 1;


Pm =  MX(:,:,21) - (MX(:,:,5)./MX(:,:,6)).^2.*(MX(:,:,20) - MX(:,:,21));
Pm =  MX(:,:,21) + (MX(:,:,5)./MX(:,:,6)).^2.*(MX(:,:,21) - MX(:,:,20));

figure();

scatter(  12*log2(MX(:,:,13)/440), log10(Pm) , 'filled'); box on; grid on;

xlabel('$12log_2(f_1/440)$','interpreter','latex'); box on; grid on;
ylabel('P_{m} [Pa] (log10)','interpreter','latex');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delays: foot and mouth rad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% t20_foot normalized by tau_f = Vf/(c_o * S_in) [HYDRODYNAMIC DELAY]

figure();
scatter( 12*log2(MX(:,:,13)/440),( co * MX(:,:,34).*MX(:,:,5)./MX(:,:,2) ), 'b', 'filled');
ylabel('$ \frac{t^{20}_f}{\tau_f} = \frac{c_o \times t^{20}_{f} \mathcal{S}_{in}}{V_f}$','interpreter','latex','Rotation',0);
xlabel('$12\times log_2(f_1/440)$','interpreter','latex'); box on; grid on;
title('Hydrodynamic delay','interpreter','latex');
ylim([0 1.3]);

%% t20_foot normalized by tau_f = Vf/(u_in * S_in) [HYDRODYNAMIC DELAY]

tau_f = MX(:,:,2)./(MX(:,:,5).*sqrt( 2/rho*(MX(:,:,20)-MX(:,:,21))) );

figure();
scatter( 12*log2(MX(:,:,13)/440), MX(:,:,34)./tau_f , 'b', 'filled');
ylabel('$ \frac{t^{20}_f}{\tau_f} = \frac{ \times t^{20}_{f} }{ \frac{V_f}{u_{in}\mathcal{S}_{in}}  }$','interpreter','latex','Rotation',0);
xlabel('$12\times log_2(f_1/440)$','interpreter','latex'); box on; grid on;
title('Hydrodynamic delay','interpreter','latex');
ylim([0 0.2]);


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
scatter(12*log2(MX(:,:,13)/440), MX(:,:,30).*MX(:,:,13), 'b','filled' ); 
xlabel('$12 \times log_2(f_1/f_{ref})$','interpreter','latex');
ylabel('$\frac{PRT_{rad}}{T_1}$','interpreter','latex','Rotation',0);
box on; grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Study of Beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% beta vs t20foot [meeh...]

figure();

scatter( MX(:,:,26), MX(:,:,34), 'b', 'filled');
grid on; box on;
xlabel('beta'); ylabel('t20foot')

%% beta vs PRTfoot [OK]
figure();
scatter( 1./MX(:,:,26), MX(:,:,29)  , 'b','filled');
xlabel('$\beta^{-1}$', 'interpreter','latex');
ylabel('PRT$_{foot}$ [s]', 'Interpreter','latex');
box on; grid on;
xlim([0 2.2e-3]);
ylim([0 4.5e-3]);

%% BETA*T1 versus NU [YES]
figure();
scatter( (MX(:,:,26)./MX(:,:,13)), (MX(:,:,27)), 'b','filled');
xlabel('beta*T1'); ylabel('nu');
grid on; box on;

%% beta vs exp(nu) [maybe]
figure();
mask = [1:10,13:22];
mask = 1:22;

scatter( (MX(:,mask,26)), (MX(:,mask,27)), 'b', 'filled');
xlabel('$\beta$', 'interpreter','latex');
ylabel('$\nu$', 'Interpreter','latex');

%% nu vs beta^(0.56) [meeh...]
figure();
mask = [1:10,13:22];
mask = 1:22;
% scatter( MX(:,mask,26), (MX(:,mask,27)), 'b', 'filled');
% scatter( log10(MX(:,mask,26))/log10(20), log(MX(:,mask,27)), 'b', 'filled');
scatter( MX(:,mask,27), MX(:,mask,26).^(0.56), 'b', 'filled');
ylabel('$\beta^{0.56}$', 'interpreter','latex');
xlabel('$\nu$', 'Interpreter','latex'); grid on; box on;





%%  Beta x T1 vs freq [OK]

mask = [1:22]; % low: sample pipes 2 and 3; high, sample pipes 11 and 12

figure();
% scatter( 12*log2(MX(:,mask,13)/440) , log10(MX(:,mask,26)./MX(:,mask,13)) ,'b','filled');ylabel('beta x T_1 (log10)');xlabel('12log_2(f_1/440)'); box on; grid on; 
scatter( 12*log2(MX(:,mask,13)/440) , (MX(:,mask,26)./MX(:,mask,13)) ,'b','filled');
ylabel('$\beta \times  T_1$','interpreter','latex');
xlabel('$12log_2(f_1/f_{ref})$','interpreter','latex'); 
box on; grid on; 
ylim([0 11]);

%% Beta x T1 vs freq BOX plot [OK]

figure();
boxplot(MX(:,:,26)./MX(:,:,13), 'Positions',12*log2(F1MEAN/440), 'Widths',1);
xlabel({'Sample pipe xticks','12log_2() fspacing'},'interpreter','latex');
ylabel('$\beta \times T_1$','interprete','latex');
ylim([0 11]);

%% nu [OK]
figure();
scatter( 12*log2(MX(:,:,13)/440) , (MX(:,:,27)) , 'b','filled'); ylabel('nu log10'); box on; grid on;

figure(); boxplot( (MX(:,:,27)) , 'Positions', 12*log2(F1MEAN/440), 'Widths', 1); ylabel('nu log10');


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



%% MaxKeyVelocity vs Anything [NO]
figure(24); clf;
for idx = [1:44]
   scatter( (MX(:,:,42)), ((MX(:,:,idx))) , 'b','filled');
   try
    xlabel('VeloMaxAbs'); ylabel([namevarsall{idx},' ']);%ylabel(num2str(idx));
   end 
   title(sprintf('Var num: %d', idx));
   box on; drawnow();
   pause();
end



%% BETA vs MaxKeyVel [NO]
figure();
scatter( MX(:,:,26), MX(:,:,42), 'b','filled');
xlabel('beta');ylabel('maxKeyVelo [m/s]');

%% NU vs MaxKeyVel [NO]
figure();
scatter( log10(1.*MX(:,:,27)), log10(1.*MX(:,:,42)), 'b','filled');
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