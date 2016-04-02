% Calculate carbon fate breakdown for varying kmH and plot bicarbonate
% transport per CO2 fixed

% Params cell with full CCM.
% 
ccm_params = CCMParams_Csome;
% 
% % set external pH
ccm_params.pH_out = 7;

% ccm_params = p;



% We are sweeping over cytoplasmic pH and calculating implied
% permeabilities of the cell membrane to bicarbonate.
% Note: we are ignoring the second pKa of carbonic acid which is at 
% about pH 10.3, so we shouldn't consider pHs too near that. Also,
% I know of no cells with an intracellular pH below 6. The only reason 
% we start the pH sweep below 6 is that the original model was using an
% implied pH of ~= 4 (calculated from the bicarbonate permeability used).
% pH = linspace(4, 8.3, 40);
% RuBisCO pH dependent data has a range of ~6 to 8.3. We extraplolate out
% to 8.9 at 8.9 the k_cat rate for carboxylation is zero, so extrapolating
% beyond that point is probably unphysical.
pH = linspace(6,8.9,100);
% pH = [7 7.33 7.5 8 8.28 8.5];
kmH = zeros(size(pH));

Hmax = 30000;   % Maximum cytoplasmic bicarbonate conc. in uM
kc_opt = 3e-5;  % optimal carboxysome permeability in kc vs jc space at pH =8
alpha = 1e-5;      % assume no conversion of cytoplasmic CO2 to bicarbonate.
ratio = 1;

% Calculate carbon fate for cells w/ carboxysomes
for i =  1:length(pH)  
    ccm_params.pH = pH(i);
    ccm_params.pH;
    ccm_params.kcC = kc_opt;
    ccm_params.kcH = ratio*kc_opt;
    ccm_params.alpha = alpha;
%     CCMParams calculates the implied kmH for us.
    kmH(i) = ccm_params.kmH_in; 
    
%     Calculate the optimal jc.
    jc_opt(i) = ccm_params.CalcOptimalJc(Hmax);
    ccm_params.jc = jc_opt(i);
    
%     Run the model
    executor = FullCCMModelExecutor(ccm_params);
    res = executor.RunAnalytical();
    Hin(i) = res.Hin_pm;
    CratewO(i) = res.CratewO_pm;
    Ccyto(i) = res.c_cyto_uM;
    Ccsome(i) = res.c_csome_uM;
    Hcsome(i) = res.h_csome_uM;
    OratewC(i) = res.OratewC_pm;
    OHrateCA(i) = res.OHrate_pm;
    Hleak(i) = res.Hleak_pm;
    Cleak(i) = res.Cleak_pm;
    Ccsomeleak(i) = res.Ccsomeleak_pm;
    
end
