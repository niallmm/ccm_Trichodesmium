% plot RuBisCO parameter pH dependence
addpath('/Users/niallmangan/GitHub/ccm/matlab')
 p = CCMParams_Csome;
rxn_Km = load('pH_Km.mat');
rxn_Vmax = load('pH_Vmax.mat');

pHvec = linspace(6, 9, 100);
p.pHmeas = 7.8;

% Whitehead parameters
% p.pHmeas = 8;
% p.Km_meas = 169;

% for color settings
c_num = 2; % number of colors on plot

p.pH_csome = p.pHmeas;
Vmaxnorm = p.kRub_pH;
Kmnorm = p.Km;
KOnorm = p.KO;

Vmaxtemp = interp1(rxn_Vmax.pH, rxn_Vmax.Vmax, p.pHmeas, 'pchip');
Kmtemp = interp1(rxn_Km.pH, rxn_Km.Km, p.pHmeas, 'pchip');

for jj = 1:length(pHvec)
    p.pH_csome = pHvec(jj);
    Vmax(jj) = p.kRub_pH;
    Km(jj) = p.Km;
    KO(jj) = p.KO;
end
    

% RuBisCO Vmax plot
figure(1) 
semilogy(rxn_Vmax.pH, rxn_Vmax.Vmax*Vmaxnorm/Vmaxtemp, 'o', 'Color', newcolor(1,c_num))
hold on
semilogy(pHvec, Vmax, '--', 'Color', newcolor(2,c_num))
xlabel('pH')
ylabel('V_{max}')

% RuBisCO Km plot
figure(2)
semilogy(rxn_Km.pH, rxn_Km.Km*Kmnorm/Kmtemp, 'o', 'Color', newcolor(1,c_num))
hold on
semilogy(pHvec, Km, '--', 'Color', newcolor(2,c_num))
xlabel('pH')
ylabel('K_m')

% RuBisCO K0 plot
figure(3) 
semilogy(rxn_Km.pH, rxn_Km.Km*KOnorm/Kmtemp, 'o', 'Color', newcolor(1,c_num))
hold on
semilogy(pHvec, KO, '--', 'Color', newcolor(2,c_num))
xlabel('pH')
ylabel('K_O')
