% sensitivity for cost

clear all
close all
clc
addpath(fileparts(pwd))
savefolder = '..\savedoutput\04242016Code\sensitivity\';
mkdir(savefolder)

p = CCMParams_Csome;
Hmax = 30000;
costpertransport = 4;


% baseline
calcCarbonfate_vs_pH_sensitivity
plotCarbonfateEnergetics2
save([savefolder, 'baseline'])
figure(2)
semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
hold on
ylabel('Energetic Cost H^+ / (CO_2 fixed)')
xlabel('Cytoplasmic pH')
drawnow
pH7_3(1) =total_cost_ccm_h(2);
pH8_3(1) =total_cost_ccm_h(4);
labels{1} = 'baseline';

% transport cost low
costpertransport = 2;
plotCarbonfateEnergetics2
save([savefolder, '2Hplus'])
figure(2)
semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
hold on
ylabel('Energetic Cost H^+ / (CO_2 fixed)')
xlabel('Cytoplasmic pH')
drawnow
pH7_3(16) =total_cost_ccm_h(2);
pH8_3(16) =total_cost_ccm_h(4);
labels{16} = '2 H+/HCO_3^-';

% transport cost high
costpertransport = 8;
plotCarbonfateEnergetics2
save([savefolder, '8Hplus'])
figure(2)
semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
hold on
ylabel('Energetic Cost H^+ / (CO_2 fixed)')
xlabel('Cytoplasmic pH')
drawnow
pH7_3(17) =total_cost_ccm_h(2);
pH8_3(17) =total_cost_ccm_h(4);
labels{17} = '8 H+/HCO_3^-';

costpertransport = 4;

% Cell membrane to CO2
clear p
p = CCMParams_Csome;
p.kmC = 0.03;
calcCarbonfate_vs_pH_sensitivity
plotCarbonfateEnergetics2
save([savefolder, 'kmC0_03'])
figure(2)
semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
drawnow
pH7_3(2) =total_cost_ccm_h(2);
pH8_3(2) =total_cost_ccm_h(4);
labels{2} = 'k_m^C 0.03 cm/s';

clear p
p = CCMParams_Csome;
p.kmC = 3;
calcCarbonfate_vs_pH_sensitivity
plotCarbonfateEnergetics2
save([savefolder, 'kmC3'])
figure(2)
semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
drawnow
pH7_3(3) =total_cost_ccm_h(2);
pH8_3(3) =total_cost_ccm_h(4);
labels{3} = 'k_m^C 3 cm/s';

% Cell membrane to H2CO3
clear p
p = CCMParams_Csome;
p.kmH_base =3e-4;
p.kcC = 1e-5;
p.kcH =1e-5;
calcCarbonfate_vs_pH_sensitivity
plotCarbonfateEnergetics2
save([savefolder, 'kmH3e_4kcC1e_5'])
figure(2)
semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
drawnow
pH7_3(4) =total_cost_ccm_h(2);
pH8_3(4) =total_cost_ccm_h(4);
labels{4} = 'k_m^{H_2CO_3} 3e-4 cm/s \newline k_c =1e-5 cm/s';

clear p
p = CCMParams_Csome;
p.kmH_base = 3e-2;
p.kcC = 1e-4;
p.kcH =1e-4;
calcCarbonfate_vs_pH_sensitivity
plotCarbonfateEnergetics2
save([savefolder, 'kmH3e_4kcC1e_4'])
semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
drawnow
pH7_3(5) =total_cost_ccm_h(2);
pH8_3(5) =total_cost_ccm_h(4);
labels{5} = 'k_m^{H_2CO_3} 3e-2 cm/s \newline k_c =1e-4 cm/s';

% % carboxysome size
% clear p
% p = CCMParams_Csome;
% p.Rc = 1e-6; % 10 nm radius
% p.kcC = 1e-3;
% p.kcH = 1e-3;
% calcCarbonfate_vs_pH_sensitivity
% plotCarbonfateEnergetics2
% save([savefolder, 'Rc10nm_kcC1e_3'])
% figure(2)
% semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
% drawnow
% pH7_3(6) =total_cost_ccm_h(2);
% pH8_3(6) =total_cost_ccm_h(4);
% labels{6} = 'R_c 10 nm \newline k_c =1e-3 cm/s';
% 
% clear p
% p = CCMParams_Csome;
% p.Rc = 10e-6; % 100 nm radius
% p.kcC = 2e-6;
% p.kcH = 2e-6;
% calcCarbonfate_vs_pH_sensitivity
% plotCarbonfateEnergetics2
% save([savefolder, 'Rc100nm2e_6'])
% semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
% drawnow
% pH7_3(7) =total_cost_ccm_h(2);
% pH8_3(7) =total_cost_ccm_h(4);
% labels{7} = 'R_c 100 nm \newline k_c =2e-6 cm/s';

%csome permeability ratio
clear p
p = CCMParams_Csome;
p.kcH = 10*p.kcC;
calcCarbonfate_vs_pH_sensitivity
plotCarbonfateEnergetics2
save([savefolder, 'ratio10'])
figure(2)
semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
drawnow
pH7_3(8) =total_cost_ccm_h(2);
pH8_3(8) =total_cost_ccm_h(4);
labels{8} = 'ratio 10 ';

clear p
p = CCMParams_Csome;
p.kcH = 100*p.kcC;
calcCarbonfate_vs_pH_sensitivity
plotCarbonfateEnergetics2
save([savefolder, 'ratio100'])
semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
pH7_3(9) =total_cost_ccm_h(2);
pH8_3(9) =total_cost_ccm_h(4);
labels{9} = 'ratio 100 ';

clear p
p = CCMParams_Csome;
p.kcH = 1000*p.kcC;
calcCarbonfate_vs_pH_sensitivity
plotCarbonfateEnergetics2
save([savefolder, 'ratio1000'])
figure(2)
semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
drawnow
pH7_3(10) =total_cost_ccm_h(2);
pH8_3(10) =total_cost_ccm_h(4);
labels{10} = 'ratio 1000 ';

%cytosolic HCO3- pool

clear p
p = CCMParams_Csome;
% p.h_cyto_exp = 5000;
Hmax =5000;
calcCarbonfate_vs_pH_sensitivity
plotCarbonfateEnergetics2
save([savefolder, 'Hcyto5mM'])
semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
pH7_3(11) =total_cost_ccm_h(2);
pH8_3(11) =total_cost_ccm_h(4);
labels{11} = 'cytosolic HCO_3^- = 5 mM ';

clear p
p = CCMParams_Csome;
% p.h_cyto_exp = 15000;
Hmax = 15000;
calcCarbonfate_vs_pH_sensitivity
plotCarbonfateEnergetics2
save([savefolder, 'Hcyto15mM'])
figure(2)
semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
drawnow
pH7_3(12) =total_cost_ccm_h(2);
pH8_3(12) =total_cost_ccm_h(4);
labels{12} =  'cytosolic HCO_3^- = 15 mM ';

clear p
p = CCMParams_Csome;
% p.h_cyto_exp = 25000;
Hmax = 25000;
calcCarbonfate_vs_pH_sensitivity
plotCarbonfateEnergetics2
save([savefolder, 'Hcyto25mM'])
figure(2)
semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
drawnow
pH7_3(13) =total_cost_ccm_h(2);
pH8_3(13) =total_cost_ccm_h(4);
labels{13} = 'cytosolic HCO_3^- = 25 mM ';

Hmax = 30000;
% carboxysome permeability away from optimum
clear p
p = CCMParams_Csome;
p.kcC = 0.02;
p.kcH = 0.02;
calcCarbonfate_vs_pH_sensitivity
plotCarbonfateEnergetics2
save([savefolder, 'kcCkcH0_02'])
figure(2)
semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
drawnow
pH7_3(14) =total_cost_ccm_h(2);
pH8_3(14) =total_cost_ccm_h(4);
labels{14} = 'k_c = 0.02 cm/s';

clear p
p = CCMParams_Csome;
p.kcC = 0.002;
p.kcH = 0.002;
calcCarbonfate_vs_pH_sensitivity
plotCarbonfateEnergetics2
save([savefolder, 'kcCkcH0_2'])
figure(2)
semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
drawnow
pH7_3(15) =total_cost_ccm_h(2);
pH8_3(15) =total_cost_ccm_h(4);
labels{15} = 'k_c = 0.002 cm/s';

figure(11)
bar([pH7_3; pH8_3]')
ylabel(
figure(12)
bar([(pH7_3-pH7_3(1))/pH7_3(1); (pH8_3-pH8_3(1))/pH8_3(1)]'*100)