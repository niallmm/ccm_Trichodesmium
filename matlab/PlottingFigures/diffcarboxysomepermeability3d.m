% script to make figure for kcH and kcC varying separately
% and jc held constant
 clear all
 close all

p = CCMParams_Csome;
% vary kcH and kcC holding the other constant at these permeabilities
% kcopt = [1e-3 3e-5 1e-6];
kcopt = 3e-5;
kcHvec = logspace(-9, 0, 100);
kcCvec = logspace(-9, 0, 100);

% calculate the CO_2 concentration in the carboxysome if there were perfect
% HCO3- permeability and perfect CO2 trapping: ie if CO2 were in
% equilibrium with HCO3- in the cytosol without concern about CO2 leakage.
% jc_opt =2.0992e-5;
% p.jc = jc_opt;
% %     Calculate the optimal jc at optimal kcC and kcH
Hmax = 30000;
p.kcC = kcopt;
p.kcH = kcopt;
jc_opt = p.CalcOptimalJc(Hmax);
p.jc = jc_opt;
    
%     Run the model
executor = FullCCMModelExecutor(p);
res = executor.RunAnalytical();

CO2max= res.h_cyto_uM*ones(size(kcHvec)) / p.Keq;

for ii = 1:length(kcHvec)
    
    for jj = 1:length(kcCvec)
        p.kcH = kcHvec(ii);
        p.kcC = kcCvec(jj);
        p.jc = jc_opt;
        
        %     Run the model
        executor = FullCCMModelExecutor(p);
        res = executor.RunAnalytical();
        Hin_c(ii,jj) = res.Hin_pm;
        CratewO_c(ii,jj) = res.CratewO_pm;
        Ccyto_c(ii,jj) = res.c_cyto_uM;
        Ccsome_c(ii,jj) = res.c_csome_uM;
        Hcsome_c(ii,jj) = res.h_csome_uM;
        OratewC_c(ii,jj) = res.OratewC_pm;
        OHrateCA_c(ii,jj) = res.OHrate_pm;
        Hleak_c(ii,jj) = res.Hleak_pm;
        Cleak_c(ii,jj) = res.Cleak_pm;
        Ccsomeleak_c(ii,jj) = res.Ccsomeleak_pm;
    end
end

figure(131);
surf(kcHvec, kcCvec, Ccsome_c);
alpha .7;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'ZScale', 'log');
% Confusingly, the major axis of the matrix is y and the minor
% is x according to surf() 
ylabel('carboxysome permeability to HCO_3^-, k_c^H');
xlabel('carboxysome permeability to CO_2, k_c^C');
zlabel('carboxysomal CO_2 concentration \muM');
