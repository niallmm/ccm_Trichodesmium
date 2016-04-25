% script to make figure for kcH and kcC varying separately
% and jc held constant
 clear all
 close all

p = CCMParams_Csome;
% vary kcH and kcC holding the other constant at these permeabilities
% kcopt = [1e-3 3e-5 1e-6];
kcopt = 3e-5;
kcHvec = logspace(-9, 0, 1e3);
kcCvec = logspace(-9, 0, 1e3);
% calculate the CO_2 concentration in the carboxysome if there were perfect
% HCO3- permeability and perfect CO2 trapping: ie if CO2 were in
% equilibrium with HCO3- in the cytosol without concern about CO2 leakage.
% jc_opt =2.0992e-5;
% p.jc = jc_opt;
    % %     Calculate the optimal jc at optimal kcC and kcH
        Hmax = 30000;
        jc_opt = p.CalcOptimalJc(Hmax);
        p.jc = jc_opt;
    
%     Run the model
executor = FullCCMModelExecutor(p);
res = executor.RunAnalytical();

CO2max= res.h_cyto_uM*ones(size(kcHvec))/p.Keq;

for jj = 1:length(kcopt)

    p.kcC = kcopt(jj);
    p.kcH = kcopt(jj);
    
    % %     Calculate the optimal jc at optimal kcC and kcH
        Hmax = 30000;
        jc_opt = p.CalcOptimalJc(Hmax);
        p.jc = jc_opt;
    
%     % hold jc at optimal value for kc = 3e-5;
%     jc_opt =2.0992e-5;
%     p.jc = jc_opt;
    
    % change permeability to CO2
    for ii = 1:length(kcCvec)
        p.kcC = kcCvec(ii);
        
        %     Run the model
        executor = FullCCMModelExecutor(p);
        res = executor.RunAnalytical();
        Hin_c(jj,ii) = res.Hin_pm;
        CratewO_c(jj,ii) = res.CratewO_pm;
        Ccyto_c(jj,ii) = res.c_cyto_uM;
        Ccsome_c(jj,ii) = res.c_csome_uM;
        Hcsome_c(jj,ii) = res.h_csome_uM;
        OratewC_c(jj,ii) = res.OratewC_pm;
        OHrateCA_c(jj,ii) = res.OHrate_pm;
        Hleak_c(jj,ii) = res.Hleak_pm;
        Cleak_c(jj,ii) = res.Cleak_pm;
        Ccsomeleak_c(jj,ii) = res.Ccsomeleak_pm;
    end
    
    p.kcC = kcopt(jj);
    % change permeability to HCO3-
    for ii = 1:length(kcHvec)
        p.kcH = kcHvec(ii);
        
        %     Run the model
        executor = FullCCMModelExecutor(p);
        res = executor.RunAnalytical();
        Hin_h(jj,ii) = res.Hin_pm;
        CratewO_h(jj,ii) = res.CratewO_pm;
        Ccyto_h(jj,ii) = res.c_cyto_uM;
        Ccsome_h(jj,ii) = res.c_csome_uM;
        Hcsome_h(jj,ii) = res.h_csome_uM;
        OratewC_h(jj,ii) = res.OratewC_pm;
        OHrateCA_h(jj,ii) = res.OHrate_pm;
        Hleak_h(jj,ii) = res.Hleak_pm;
        Cleak_h(jj,ii) = res.Cleak_pm;
        Ccsomeleak_h(jj,ii) = res.Ccsomeleak_pm;
    end
end
for jj = 1:length(kcopt)
% tag = {'--','-', ':'};
tag = {'-'};
figure(131)
loglog(kcCvec, Ccsome_c(jj,:), [tag{jj} 'r'])
hold on
loglog(kcHvec, Ccsome_h(jj,:), [tag{jj} 'b'])
xlabel('carboxysome permeability')
ylabel('CO_2 concentration in carboxysome')

end
plot(kcCvec, CO2max, '--k')
    legend(    'k_c^C = 3\times10^{-5} cm/s, k_c^H varying',...
    'k_c^H = 3\times10^{-5} cm/s, k_c^C varying',...
        'CO_2 in equilibrium with HCO_3^-', 'Location', 'southeast')
legend('boxoff')


costpertransport = 4;
% Total cost of the system in H+/fixation for each case.
total_cost_ccm_c = totalProtonCost_CCM(Hin_c, CratewO_c, OratewC_c, ...
                    OHrateCA_c, costpertransport);
total_cost_ccm_h = totalProtonCost_CCM(Hin_h, CratewO_h, OratewC_h, ...
                    OHrateCA_h, costpertransport);
                

clear p
p = CCMParams_Csome;
p.kcC = 1e-16;
p.kcH = 1e15;
% %     Calculate the optimal jc at optimal kcC and kcH
Hmax = 30000;
jc_opt = p.CalcOptimalJc(Hmax);
p.jc = jc_opt;

executor = FullCCMModelExecutor(p);
res = executor.RunAnalytical();
Hin_min = res.Hin_pm;
CratewO_min = res.CratewO_pm;
Ccyto_min = res.c_cyto_uM;
Ccsome_min = res.c_csome_uM;
Hcsome_min = res.h_csome_uM;
OratewC_min = res.OratewC_pm;
OHrateCA_min = res.OHrate_alt_pm;
Hleak_min = res.Hleak_pm;
Cleak_min = res.Cleak_pm;
Ccsomeleak_min = res.Ccsomeleak_pm;
% have to set cost of OHrateCA to zero, because we cant check leakage of
% CO2 out of carboxysome anymore.
total_cost_ccm_min = totalProtonCost_CCM(Hin_min, CratewO_min, OratewC_min, ...
                   OHrateCA_min, costpertransport);
                
        
tag = {'-'};
figure(132)
loglog(kcHvec, total_cost_ccm_h, [tag{1} 'b'])
hold on
plot(kcHvec, total_cost_ccm_min*ones(size(kcHvec)), '--k')
xlabel('carboxysome permeability')
ylabel('Cost of CCM [H^+/(CO_2 fixation)]')


    legend(    'k_c^C = 3\times10^{-5} cm/s, k_c^H varying',...
        'CO_2 in equilibrium with HCO_3^-', 'Location', 'southeast')
legend('boxoff')