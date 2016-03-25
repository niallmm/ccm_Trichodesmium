% plot RuBisCO parameter pH dependence
addpath('/Users/niallmangan/GitHub/ccm/matlab')
 p = CCMParams_Csome;
rxn_Km = load('pH_Km.mat');
rxn_Vmax = load('pH_Vmax.mat');

pHvec = linspace(6, 9, 100);

% for color settings
c_num = 2; % number of colors on plot


for jj = 1:length(pHvec)
    p.pH = pHvec(jj);
    Vmax(jj) = p.kRub_pH;
    Km(jj) = p.Km;
    KO(jj) = p.KO;
end
    

% RuBisCO Vmax plot
figure(1) 
semilogy(rxn_Vmax.pH, rxn_Vmax.Vmax*p.kRub/2.1702, 'o', 'Color', newcolor(1,c_num))
hold on
semilogy(pHvec, Vmax, '--', 'Color', newcolor(2,c_num))

% RuBisCO Km plot
figure(2)
semilogy(rxn_Km.pH, rxn_Km.Km*1e3, 'o', 'Color', newcolor(1,c_num))
hold on
semilogy(pHvec, Km, '--', 'Color', newcolor(2,c_num))

% RuBisCO K0 plot
figure(3) 
semilogy(rxn_Km.pH, rxn_Km.Km*p.KO7_8/p.Km7_8*1e3, 'o', 'Color', newcolor(1,c_num))
hold on
semilogy(pHvec, KO, '--', 'Color', newcolor(2,c_num))
