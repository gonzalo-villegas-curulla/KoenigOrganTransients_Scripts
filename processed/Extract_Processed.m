clc, clear;

try 
    run ../Analysis_PlotProcessedData_Loads.m
catch
    run Analysis_PlotProcessedData_Loads.m
end

rho = 1.2; % Air density [kg/m^3]
co  = 340; % Acoustic propagation speed [m/s]
co2 = co^2; 
P0  = 820; % Pallet box pressure [Pa = kg m^-1 s^-2]

Vgrv = median(MX(:,:,11),1,'omitnan');
Vf   = median(MX(:,:,2),1, 'omitnan');

Spall_Slot    = median(MX(:,:,3),1,'omitnan')*0.1298; % Perforated rectangles on the plate, with the same width as the groove channel
Spall_Lateral = PalletValveStrokeArea(maskpipes); % At maximum aperture of valve, adding areas of a rectangle and two triangles
Spall_eff_meas = median(MX(:,:,6),1,'omitnan').*sqrt(median(MX(:,:,21),1,'omitnan'))./sqrt( median(MX(:,:,19),1,'omitnan')-median(MX(:,:,20),1,'omitnan'));


MX = squeeze(median(MX,1,'omitnan'))';

pipelist = [3,4,5,6,7,9,10,11,13,15,17,19,24,25,27,29,32,34,37,39,41,44];

data_proc = {};

data_proc.Lp      = MX(1,:);
data_proc.Vf      = MX(2,:);
data_proc.Vgrv    = MX(11,:);

data_proc.PWidth  = MX(3,:);


data_proc.Spall_geom = Spall_Lateral';
data_proc.Sslot_geom = Spall_Slot;
data_proc.Sgrv_geom = MX(3,:)*0.05;
data_proc.Sin_geom     = MX(5,:);
data_proc.Sjet_geom    = MX(6,:);
data_proc.Spall_eff = Spall_eff_meas;


data_proc.F1      = MX(13,:);
data_proc.Ppall   = MX(19,:);
data_proc.Pgrv    = MX(20,:);
data_proc.Pf      = MX(21,:);
data_proc.Prad    = MX(22,:);

data_proc.Beta    = MX(26,:);

data_proc.t10grv  = MX(51,:);
data_proc.PRTgrv  = MX(52,:);
data_proc.t10ft   = MX(55,:);
data_proc.PRT10ft = MX(56,:);
data_proc.t10rad  = MX(57,:);
data_proc.PRT10rad  = MX(58,:);

% -------------------------------

Sgrv_eff = MX(3,:)'*0.05;

% % % % Spall_max = PalletValveStrokeArea(maskpipes);
% % % %     Spall_max = Spall_max(:); Spall_max = Spall_max(:);
% % % % Sin           = MX(5,:);
% % % %     Sin = median(Sin, 1, 'omitnan'); Sin = Sin(:);
Sj            = 1*MX(6,:);
    Sj = median(Sj, 1, 'omitnan'); Sj = Sj(:);


        % % % % Spall_max_eff = 1.*Spall_max;
Spall_eff = 10.^(  -0.03732*12*log2(F1MEAN/440) - 4.6202 ); % fitted to computed 
Sin_eff = Sj'.*sqrt( MX(21,:)./(MX(19,:)-MX(21,:)) );
    Sin_eff = Sin_eff';
Sj_eff = 1.*Sj;

Vf_over_Vgrv = Vf./Vgrv;

one_over_Amax = Spall_eff*co2./Vgrv'*sqrt(rho/P0);
one_over_B    = Sin_eff*co2./Vgrv'*sqrt(rho/P0);
one_over_C    = Sin_eff*co2./Vf'*sqrt(rho/P0);
one_over_D    = Sj_eff*co2./Vf'*sqrt(rho/P0);
sigMa         = Spall_eff./Sgrv_eff;
%----------------------------------------
data_proc.one_over_Amax = one_over_Amax ;
data_proc.one_over_B = one_over_B;
data_proc.one_over_C = one_over_C;
data_proc.one_over_D = one_over_D;
data_proc.sigma = sigMa;


save('data_proc.mat','data_proc');
save('./NumericalModel/data_proc.mat','data_proc');
