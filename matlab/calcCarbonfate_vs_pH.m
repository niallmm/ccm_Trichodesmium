% Calculate carbon fate breakdown for varying kmH and plot bicarbonate
% transport per CO2 fixed

% Params cell with full CCM.

ccm_params = CCMParams_Csome;

% set external pH
ccm_params.pH_out = 7;



% We are sweeping over cytoplasmic pH and calculating implied
% permeabilities of the cell membrane to bicarbonate.
% Note: we are ignoring the second pKa of carbonic acid which is at 
% about pH 10.3, so we shouldn't consider pHs too near that. Also,
% I know of no cells with an intracellular pH below 6. The only reason 
% we start the pH sweep below 6 is that the original model was using an
% implied pH of ~= 4 (calculated from the bicarbonate permeability used).
% pH = linspace(4, 8.3, 40);
pH = linspace(6,8.3,40);
kmH = zeros(30);

Hmax = 30000;   % Maximum cytoplasmic bicarbonate conc. in uM
kc_opt = 1e-4;  % optimal carboxysome permeability from previous work.
alpha = 0;      % assume no conversion of cytoplasmic CO2 to bicarbonate.

% Calculate carbon fate for cells w/ carboxysomes
for i =  1:length(pH)  
    ccm_params.pH = pH(i);
    ccm_params.kcC = kc_opt;
    ccm_params.kcH = kc_opt;
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
    
end
