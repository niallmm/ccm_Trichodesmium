% Calculate carbon fate breakdown for varying kmH and plot bicarbonate
% transport per CO2 fixed


% Params for 2 cases.
% 1: cell has full CCM.
% 2: cell has just RuBisCO.
ccm_params = CCMParams_Csome;
ccm_params_just_rbc = CCMParams_NoCsome;

% ccm_params.kRub = 11.6;   % rxns/s maximum reaction rate at single active site
% ccm_params.Km_8 = 340;    % half max reaction rate of RuBisCO, uM
% ccm_params.S_sat = 43;    % specificity ratio
% ccm_params.KO_8 = 972;    % uM

% As a thought experiment, use RuBisCO parameters from Goldiera Sulfuraria
% which has a much higher specificity to CO2 and much lower KM than the
% cyanobacterial RuBisCO. Data from Savir et al., 2010.
ccm_params_just_rbc.kRub = 1.2;
ccm_params_just_rbc.Km7_8 = 3.3; 
ccm_params_just_rbc.S_sat = 166;
ccm_params_just_rbc.KO7_8 = 374;
% 
% ccm_params_just_rbc.kRub = 11.6;
% ccm_params_just_rbc.Km7_8 = 340; 
% ccm_params_just_rbc.S_sat = 43;
% ccm_params_just_rbc.KO7_8 = 972;

ccm_params_just_rbc.O = 300;
ccm_params.O = 300;
% 
% ccm_params_just_rbc.I_in = 0;
% ccm_params_just_rbc.I_out = 0;
% ccm_params.I_in = 0;
% ccm_params.I_out = 0;


 % set exeternal pH
ccm_params.pH_out = 7;
ccm_params_just_rbc.pH_out = 7;

Hmax = 30000;   % Maximum cytoplasmic bicarbonate conc. in uM
kc_opt = 1e-4;  % optimal carboxysome permeability from previous work.
ccm_params.kcH = kc_opt;
ccm_params.kcC = kc_opt;
alpha = 0;      % assume no conversion of cytoplasmic CO2 to bicarbonate.

ccm_params.pH = 8.3; % internal pH
ccm_params_just_rbc.pH = 8.3;


CO2extv = logspace(-2, 5, 50);


% Calculate carbon fate for cells w/ carboxysomes
for i =  1:length(CO2extv)  
  
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




% Calculate carbon fate for cells w/o just RuBisCO
% as our approximation assumes flux in dominated by active transport.

    % set permeability kmH
    ccm_params_just_rbc.jc = 0;     % no active transport
    ccm_params_just_rbc.NCA = 1e-7; % no carbonic anhydrase
    
for i =  1:length(CO2extv)
    
     ccm_params_just_rbc.Cout = CO2extv(i);
    
    % note: kc meaningless in this case as no cbsome
    % note: using > 0 NCA for now to avoid some stupid error in underlying
    % code.
    
    executor = NoCsomeModelExecutor(ccm_params_just_rbc);
    res = executor.RunAnalytical();
    Hin_just_rbc(i) = res.Hin_pm;
    CratewO_just_rbc(i) = res.CratewO_pm;
    Ccyto_just_rbc(i) = res.c_cyto_uM;
    OratewC_just_rbc(i) = res.OratewC_pm;
end

