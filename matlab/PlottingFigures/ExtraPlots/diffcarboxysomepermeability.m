% script to make figure for kcH and kcC varying separately
% and jc held constant
% this script calculates the CO2 concentration for single jc and more kcC
% an kcH values than in paper.
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
tag = {'--','-', ':'};
figure(131)
loglog(kcCvec, Ccsome_c(jj,:), [tag{jj} 'r'])
hold on
loglog(kcHvec, Ccsome_h(jj,:), [tag{jj} 'b'])
xlabel('carboxysome permeability')
ylabel('CO_2 concentration in carboxysome')

end
plot(kcCvec, CO2max, '--k')
legend('k_c^C = 10^{-3} cm/s, k_c^H varying',...
    'k_c^H = 10^{-3} cm/s, k_c^C varying',...
    'k_c^C = 3\times10^{-5} cm/s, k_c^H varying',...
    'k_c^H = 3\times10^{-5} cm/s, k_c^C varying',...
    'k_c^C = 10^{-6} cm/s, k_c^H varying',...    
    'k_c^H = 10^{-6} cm/s, k_c^C varying',...
    'CO_2 in equilibrium with HCO_3^-', 'Location', 'southeast')
legend('boxoff')


costpertransport = 4;
% Total cost of the system in H+/fixation for each case.
total_cost_ccm_c = totalProtonCost_CCM(Hin_c, CratewO_c, OratewC_c, ...
                    OHrateCA_c, costpertransport);
total_cost_ccm_h = totalProtonCost_CCM(Hin_h, CratewO_h, OratewC_h, ...
                    OHrateCA_h, costpertransport);
                


for jj = 1:length(kcopt)
tag = {'--','-', ':'};
figure(132)
loglog(kcCvec, total_cost_ccm_c(jj,:), [tag{jj} 'r'])
hold on
loglog(kcHvec, total_cost_ccm_h(jj,:), [tag{jj} 'b'])
xlabel('carboxysome permeability')
ylabel('Cost of CCM [H^+/(CO_2 fixation)]')

end
legend('k_c^C = 10^{-3} cm/s, k_c^H varying',...
    'k_c^H = 10^{-3} cm/s, k_c^C varying',...
    'k_c^C = 3\times10^{-5} cm/s, k_c^H varying',...
    'k_c^H = 3\times10^{-5} cm/s, k_c^C varying',...
    'k_c^C = 10^{-6} cm/s, k_c^H varying',...    
    'k_c^H = 10^{-6} cm/s, k_c^C varying', 'Location', 'southeast')
legend('boxoff')