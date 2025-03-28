clc, clear;

try 
    run ../Analysis_PlotProcessedData_Loads.m
catch
    run Analysis_PlotProcessedData_Loads.m
end

data_proc = {};

% =============================
rho = 1.2; % Air density [kg/m^3]
co  = 340; % Acoustic propagation speed [m/s]
co2 = co^2; 
P0  = 820; % Pallet box pressure [Pa = kg m^-1 s^-2]    
% ==========================
Vgrv = median(MX(:,:,11),1,'omitnan');
Vf   = median(MX(:,:,2),1, 'omitnan');

Spall_Slot     = median(MX(:,:,3),1,'omitnan')*0.1298; % Perforated rectangles on the plate, with the same width as the groove channel
Spall_Lateral  = PalletValveStrokeArea(maskpipes); % At maximum aperture of valve, adding areas of a rectangle and two triangles
Spall_eff_meas = median(MX(:,:,6),1,'omitnan').*sqrt(median(MX(:,:,21),1,'omitnan'))./sqrt( median(MX(:,:,19),1,'omitnan')-median(MX(:,:,20),1,'omitnan'));
% ==================================

% For meas-simulation comparison (before flattening MX)
data_proc.Ppall_mean = mean(MX(:,:,19),1,'omitnan'); % <P0>
data_proc.Ppall_std  = std(MX(:,:,19), [], 1, 'omitnan');
data_proc.Pgrv_mean  = mean(MX(:,:,20), 1, 'omitnan');
data_proc.Pgrv_std   = std(MX(:,:,20), [], 1, 'omitnan');
data_proc.Pf_mean    = mean(MX(:,:,21), 1, 'omitnan');
data_proc.Pf_std     = std(MX(:,:,21), [], 1, 'omitnan');

data_proc.PRTgrv_mean = mean(MX(:,:,52), 1, 'omitnan');
data_proc.PRTgrv_std  = std(MX(:,:,52), [], 1, 'omitnan');
data_proc.PRTf_mean   = mean(MX(:,:,56), 1, 'omitnan');
data_proc.PRTf_std    = std(MX(:,:,56), [], 1, 'omitnan');

%==================================

MXbuf = MX;
MX    = squeeze(median(MX,1,'omitnan'))';

pipelist = [3,4,5,6,7,9,10,11,13,15,17,19,24,25,27,29,32,34,37,39,41,44];


data_proc.Lp      = MX(1,:);
data_proc.Vf      = MX(2,:);
data_proc.Vgrv    = MX(11,:);
data_proc.PWidth  = MX(3,:);

%==================================

data_proc.Spall_geom = Spall_Lateral';
data_proc.Sslot_geom = Spall_Slot;
data_proc.Spall_eff  = Spall_eff_meas;

data_proc.Sgrv_geom  = MX(3,:)*0.05;

data_proc.Sin_geom   = MX(5,:);

data_proc.Sjet_geom  = MX(6,:);



data_proc.F1       = MX(13,:);
data_proc.Ppall    = MX(19,:);
data_proc.Pgrv     = MX(20,:);
data_proc.Pf       = MX(21,:);
data_proc.Prad     = MX(22,:);

data_proc.Beta     = MX(26,:);

data_proc.t10grv   = MX(51,:);
data_proc.PRTgrv   = MX(52,:);

data_proc.t10ft    = MX(55,:);
data_proc.PRTft    = MX(56,:);

data_proc.t10rad   = MX(57,:);
data_proc.PRTrad   = MX(58,:);

% -------------------------------

Sgrv_eff = MX(3,:)'*0.05;

% % % % Spall_max = PalletValveStrokeArea(maskpipes);
% % % %     Spall_max = Spall_max(:); Spall_max = Spall_max(:);
% % % % Sin           = MX(5,:);
% % % %     Sin = median(Sin, 1, 'omitnan'); Sin = Sin(:);
Sj            = 1*MX(6,:);
    Sj = median(Sj, 1, 'omitnan'); Sj = Sj(:);


Spall_eff = 10.^(  -0.03732*12*log2(F1MEAN/440) - 4.6202 ); % fitted to computed 
% Sin_eff = Sj'.*sqrt( MX(21,:)./(MX(19,:)-MX(21,:)) ); % NOO! Noo! No!
Sin_eff = Sj'.*sqrt( MX(21,:)./(MX(20,:)-MX(21,:)) );
    Sin_eff = Sin_eff';
Sj_eff = 1.*Sj;

Vf_over_Vgrv = Vf./Vgrv;

%----------------------------------------

data_proc.Amax  = Vgrv(:)./Spall_eff./co2.*sqrt(median(MXbuf(:,:,19),1,'omitnan')'/rho);
data_proc.B     = Vgrv(:)./Sin_eff./co2.*sqrt(median(MXbuf(:,:,19),1,'omitnan')'/rho);
data_proc.C     = Vf(:)./Sin_eff./co2.*sqrt(median(MXbuf(:,:,19),1,'omitnan')'/rho);
data_proc.D     = Vf(:)./Sj_eff./co2.*sqrt(median(MXbuf(:,:,19),1,'omitnan')'/rho);
data_proc.sigma = Spall_eff./Sgrv_eff; % Sigma max







save('data_proc.mat','data_proc');
save('./NumericalModel/data_proc.mat','data_proc');
