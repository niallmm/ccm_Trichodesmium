% Calculate carbon fate breakdown for varying kmH and plot bicarbonate
% transport per CO2 fixed

% Params for 3 cases.
% 1: cell has full CCM.
% 2: cell has active transport and CA, no carboxysome.
% 3: cell has just RuBisCO.
ccm_params = CCMParams_Csome;
ccm_params_scaffold = CCMParams_Csome;
ccm_params_cell = CCMParams_NoCsome;
ccm_params_just_rbc = CCMParams_NoCsome;

% As a thought experiment, use RuBisCO parameters from Goldiera Sulfuraria
% which has a much higher specificity to CO2 and much lower KM than the
% cyanobacterial RuBisCO. Data from Savir et al., 2010.
ccm_params_just_rbc.kRub = 1.2;
ccm_params_just_rbc.Km = 3.3; 
ccm_params_just_rbc.S_sat = 166;
ccm_params_just_rbc.KO = 374;

% We are sweeping over cytoplasmic pH and calculating implied
% permeabilities of the cell membrane to bicarbonate.
% Note: we are ignoring the second pKa of carbonic acid which is at 
% about pH 10.3, so we shouldn't consider pHs too near that. Also,
% I know of no cells with an intracellular pH below 6. The only reason 
% we start the pH sweep below 6 is that the original model was using an
% implied pH of ~= 4 (calculated from the bicarbonate permeability used).
Houtv = logspace(1, 4, 40);
Coutv = 0.01*Houtv;
kmH = zeros(30);

Hmax = 30000;   % Maximum cytoplasmic bicarbonate conc. in uM
kc_opt = 1e-4;  % optimal carboxysome permeability from previous work.
alpha = 0;      % assume no conversion of cytoplasmic CO2 to bicarbonate.

ccm_params.pH = 8;
ccm_params.k = kc_opt;
% Calculate carbon fate for cells w/ carboxysomes
for i =  1:length(Houtv)  
    ccm_params.Hout = Houtv(i);
    ccm_params.Cout = Coutv(i);
    
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
end

% Calculate carbon fate for cells w/ carboxysomes but no shell
% TODO: warning that CA not saturated in this case?? does this make sense?

    ccm_params_scaffold.pH = 8;
    ccm_params_scaffold.k = 1.0;
for i =  1:length(Houtv)  
    ccm_params_scaffold.Hout = Houtv(i);
    ccm_params_scaffold.Cout = Coutv(i);
    % Calculate the optimal jc.
    jc_opt_scaffold(i) = ccm_params_scaffold.CalcOptimalJc(Hmax);
    ccm_params_scaffold.jc = jc_opt(i);
    
    % Run the model
    executor = FullCCMModelExecutor(ccm_params_scaffold);
    res = executor.RunAnalytical();
    Hin_scaffold(i) = res.Hin_pm;
    CratewO_scaffold(i) = res.CratewO_pm;
    Ccyto_scaffold(i) = res.c_cyto_uM;
    Ccsome_scaffold(i) = res.c_csome_uM;
    Hcsome_scaffold(i) = res.h_csome_uM;
end

    ccm_params_cell.pH = 8;

% Calculate carbon fate for cells w/o carboxysomes
for i =  1:length(Houtv)
    ccm_params_cell.Hout = Houtv(i);
    ccm_params_cell.Cout = Coutv(i);
    % CCMParams calculates the implied kmH for us.
    jc_opt_cell(i) = ccm_params_cell.CalcOptimalJc(Hmax);
    ccm_params_cell.jc = jc_opt_cell(i);
    
    % note: k (ci uptake rate into cbsome) meaningless in this case as 
    % there is no no cbsome
    
    % Run the model
    executor = NoCsomeModelExecutor(ccm_params_cell);
    res = executor.RunAnalytical();
    Hin_cy(i) = res.Hin_pm;
    CratewO_cy(i) = res.CratewO_pm;
    Ccyto_cy(i) = res.c_cyto_uM;
    OratewC_cy(i) = res.OratewC_pm;
end

% Calculate carbon fate for cells w/o just RuBisCO
% TODO(flamholz): may need to account for differt pH inside and outside
% as our approximation assumes flux in dominated by active transport.

    % set permeability kmH
    ccm_params_just_rbc.pH = 8;
    ccm_params_just_rbc.jc = 0;  % no active transport
    ccm_params_just_rbc.NCA = 1e-3; % no carbonic anhydrase
    
for i =  1:length(Houtv)
    
    ccm_params_just_rbc.Hout = Houtv(i);
    ccm_params_just_rbc.Cout = Coutv(i);
    
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

