% Calculate carbon fate breakdown for varying kmH and plot bicarbonate
% transport per CO2 fixed

% Params for 3 cases.
% 1: cell has full CCM.
% 2: cell has active transport and CA, no carboxysome.
% 3: cell has just RuBisCO.
ccm_params = CCMParams_Csome;
ccm_params_cell = CCMParams_NoCsome;
ccm_params_just_rbc = CCMParams_NoCsome;

% We are sweeping over cytoplasmic pH and calculating implied
% permeabilities of the cell membrane to bicarbonate.
pH = linspace(2, 9.5, 30);
kmH = zeros(30);

Hmax = 30000;   % Maximum cytoplasmic bicarbonate conc. in uM
kc_opt = 1e-4;  % optimal carboxysome permeability from previous work.
alpha = 0;      % assume no conversion of cytoplasmic CO2 to bicarbonate.

% Calculate carbon fate for cells w/ carboxysomes
for i =  1:length(pH)  
    ccm_params.pH = pH(i);
    ccm_params.k = kc_opt;
    % CCMParams calculates the implied kmH for us.
    kmH(i) = ccm_params.kmH; 
    
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
end

% Calculate carbon fate for cells w/o carboxysomes
for i =  1:length(pH)
    ccm_params_cell.pH = pH(i);
    
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
end

% Calculate carbon fate for cells w/o just RuBisCO
% TODO(flamholz): may need to account for differt pH inside and outside
% as our approximation assumes flux in dominated by active transport.
for i =  1:length(pH)
    % set permeability kmH
    ccm_params_just_rbc.pH = pH(i);
    ccm_params_just_rbc.jc = 0;  % no active transport
    ccm_params_just_rbc.NCA = 1e-3; % no carbonic anhydrase
    % note: kc meaningless in this case as no cbsome
    % note: using > 0 NCA for now to avoid some stupid error in underlying
    % code.
    
    executor = NoCsomeModelExecutor(ccm_params_just_rbc);
    res = executor.RunAnalytical();
    Hin_just_rbc(i) = res.Hin_pm;
    CratewO_just_rbc(i) = res.CratewO_pm;
    Ccyto_just_rbc(i) = res.c_cyto_uM;
end

