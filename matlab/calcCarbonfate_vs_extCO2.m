% Calculate carbon fate breakdown for varying kmH and plot bicarbonate
% transport per CO2 fixed

% Params for 3 cases.
% 1: cell has full CCM.
% 2: cell has no CCM and C3 plant RuBisCO
% 3: cell has no CCM and highly specific algal RuBisCO
ccm_params = CCMParams_Csome;
ccm_params_high_perm = CCMParams_Csome;
ccm_params_low_pH = CCMParams_Csome;
params_just_c3_rbc = CCMParams_NoCsome;
params_just_specific_rbc = CCMParams_NoCsome;

% As a thought experiment, use RuBisCO parameters from Tobacco.
% Tobacco RuBisCO is intermediate between the cyanobacterial enzyme (high
% rate, low specificity) and the Goldiera Sulfuraria enzyme (below, low
% rate and high specificity). Data from Savir et al., 2010.
params_just_c3_rbc.kRub = 3.4;      % rxns/s maximum reaction rate at single active site
params_just_c3_rbc.Km7_8 = 10.7;    % half max concentration for CO2, uM
params_just_c3_rbc.S_sat = 82;      % specificity ratio
params_just_c3_rbc.KO7_8 = 295;     % half max concentration for O2, uM

% As a thought experiment, use RuBisCO parameters from Goldiera Sulfuraria
% which has a much higher specificity to CO2 and much lower KM than the
% cyanobacterial RuBisCO. Data from Savir et al., 2010.
params_just_specific_rbc.kRub = 1.2;
params_just_specific_rbc.Km7_8 = 3.3; 
params_just_specific_rbc.S_sat = 166;
params_just_specific_rbc.KO7_8 = 374;

% set exeternal pH
ccm_params.pH_out = 7;
ccm_params_high_perm.pH_out = 7;
ccm_params_low_pH.pH_out = 7;
params_just_c3_rbc.pH_out = 7;
params_just_specific_rbc.pH_out = 7;

ccm_params.pH = 8.3; % internal pH
ccm_params_high_perm.pH = 8.3;
ccm_params_low_pH.pH = 7.3;  % cytosolic pH measured in dark
params_just_c3_rbc.pH = 8.3;
params_just_specific_rbc.pH = 8.3;

Hmax = 30000;   % Maximum cytoplasmic bicarbonate conc. in uM
kc_opt = 3e-5;  % optimal carboxysome permeability from jc kc plot at pH=8.
ccm_params.kcH = kc_opt;
ccm_params.kcC = kc_opt;
ccm_params_low_pH.kcH = kc_opt;
ccm_params_low_pH.kcC = kc_opt;
% highest feasible csome permeability from phase diagram
ccm_params_high_perm.kcH = 1e-2;  
ccm_params_high_perm.kcC = 1e-2;
alpha = 0;      % assume no conversion of cytoplasmic CO2 to bicarbonate.


CO2extv = logspace(-1, 2, 1e3);

disp('full ccm');
% Calculate carbon fate for cells w/ carboxysomes
for i = 1:length(CO2extv)  
    ccm_params.Cout = CO2extv(i);
    
    % Calculate the optimal jc.
    jc_opt(i) = ccm_params.CalcOptimalJc(Hmax);
    ccm_params.jc = jc_opt(i);
    
    % Run the model
    executor = FullCCMModelExecutor(ccm_params);
    res = executor.RunAnalytical();
    Hin(i) = res.Hin_pm;
    CratewO(i) = res.CratewO_pm;
    Ccyto(i) = res.c_cyto_uM;
    Ccsome(i) = res.c_csome_uM;
    Hcsome(i) = res.h_csome_uM;
    OratewC(i) = res.OratewC_pm;
    OHrateCA(i) = res.OHrate_pm;
    Houtv(i) = ccm_params.Hout;
end

disp('full ccm, high permeability csome');
% Calculate carbon fate for cells w/ carboxysomes
for i = 1:length(CO2extv)  
    ccm_params_high_perm.Cout = CO2extv(i);
    
    % Calculate the optimal jc.
    jc_opt(i) = ccm_params_high_perm.CalcOptimalJc(Hmax);
    ccm_params_high_perm.jc = jc_opt(i);
    
    % Run the model
    executor = FullCCMModelExecutor(ccm_params_high_perm);
    res = executor.RunAnalytical();
    Hin_high_perm(i) = res.Hin_pm;
    CratewO_high_perm(i) = res.CratewO_pm;
    Ccyto_high_perm(i) = res.c_cyto_uM;
    Ccsome_high_perm(i) = res.c_csome_uM;
    Hcsome_high_perm(i) = res.h_csome_uM;
    OratewC_high_perm(i) = res.OratewC_pm;
    OHrateCA_high_perm(i) = res.OHrate_pm;
    Houtv_high_perm(i) = ccm_params_high_perm.Hout;
end

disp('full ccm, low cytosolic pH');
% Calculate carbon fate for cells w/ carboxysomes
for i = 1:length(CO2extv)  
    ccm_params_low_pH.Cout = CO2extv(i);
    
    % Calculate the optimal jc.
    jc_opt(i) = ccm_params_low_pH.CalcOptimalJc(Hmax);
    ccm_params_low_pH.jc = jc_opt(i);
    
    % Run the model
    executor = FullCCMModelExecutor(ccm_params_low_pH);
    res = executor.RunAnalytical();
    Hin_low_pH(i) = res.Hin_pm;
    CratewO_low_pH(i) = res.CratewO_pm;
    Ccyto_low_pH(i) = res.c_cyto_uM;
    Ccsome_low_pH(i) = res.c_csome_uM;
    Hcsome_low_pH(i) = res.h_csome_uM;
    OratewC_low_pH(i) = res.OratewC_pm;
    OHrateCA_low_pH(i) = res.OHrate_pm;
    Houtv_low_pH(i) = ccm_params_low_pH.Hout;
end

% Calculate carbon fate for cells w/o just RuBisCO
% as our approximation assumes flux in dominated by active transport.

% set permeability kmH
params_just_c3_rbc.jc = 0;     % no active transport
params_just_c3_rbc.NCA = 1e-7; % no carbonic anhydrase
disp('just c3 rubisco');

for i = 1:length(CO2extv)
    params_just_c3_rbc.Cout = CO2extv(i);
    
    % note: kc meaningless in this case as no cbsome
    % note: using > 0 NCA for now to avoid some stupid error in underlying
    % code.
    
    executor = NoCsomeModelExecutor(params_just_c3_rbc);
    res = executor.RunAnalytical();
    Hin_just_c3_rbc(i) = res.Hin_pm;
    CratewO_just_c3_rbc(i) = res.CratewO_pm;
    Ccyto_just_c3_rbc(i) = res.c_cyto_uM;
    OratewC_just_c3_rbc(i) = res.OratewC_pm;
end


% set permeability kmH
params_just_specific_rbc.jc = 0;     % no active transport
params_just_specific_rbc.NCA = 1e-7; % no carbonic anhydrase
disp('just specific rubisco');

for i =  1:length(CO2extv)
    params_just_specific_rbc.Cout = CO2extv(i);
    
    % note: kc meaningless in this case as no cbsome
    % note: using > 0 NCA for now to avoid some stupid error in underlying
    % code.
    
    executor = NoCsomeModelExecutor(params_just_specific_rbc);
    res = executor.RunAnalytical();
    Hin_just_specific_rbc(i) = res.Hin_pm;
    CratewO_just_specific_rbc(i) = res.CratewO_pm;
    Ccyto_just_specific_rbc(i) = res.c_cyto_uM;
    OratewC_just_specific_rbc(i) = res.OratewC_pm;
end

