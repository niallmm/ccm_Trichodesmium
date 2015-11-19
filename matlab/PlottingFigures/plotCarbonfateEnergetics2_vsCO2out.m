% Run calcCarbonfate_vs_extCO2.m first to populate the workspace.
all_ones = ones(size(CO2extv));

% Total cost of the system in H+/fixation for each case.
total_cost_ccm_h = totalProtonCost_CCM(Hin, CratewO, OratewC, OHrateCA);
total_cost_low_pH_ccm_h = totalProtonCost_CCM(Hin_low_pH, ...
    CratewO_low_pH, OratewC_low_pH, OHrateCA_low_pH);
total_cost_high_perm_ccm_h = totalProtonCost_CCM(Hin_high_perm, ...
    CratewO_high_perm, OratewC_high_perm, OHrateCA_high_perm);

total_cost_just_c3_rbc_h = totalProtonCost_NoCCM(CratewO_just_c3_rbc, ...
    OratewC_just_c3_rbc);
total_cost_just_specific_rbc_h = totalProtonCost_NoCCM(...
    CratewO_just_specific_rbc, OratewC_just_specific_rbc); 

min_total_ccm_cost = bsxfun(@min, total_cost_ccm_h, total_cost_low_pH_ccm_h); 
min_total_ccm_cost = bsxfun(@min, min_total_ccm_cost, total_cost_high_perm_ccm_h);
max_total_ccm_cost = bsxfun(@max, total_cost_ccm_h, total_cost_low_pH_ccm_h); 
max_total_ccm_cost = bsxfun(@max, max_total_ccm_cost, total_cost_high_perm_ccm_h); 

% Absolute rates as a function of CO2. 
figure(333);

loglog(CO2extv, CratewO, 'k');
hold on
loglog(CO2extv, OratewC, '--k');
loglog(CO2extv, CratewO_just_c3_rbc, 'g');
loglog(CO2extv, OratewC_just_c3_rbc, '--g');
loglog(CO2extv, CratewO_just_specific_rbc, 'b');
loglog(CO2extv, OratewC_just_specific_rbc, '--b');
%axis([1e-2 1e5 10 4e5]);
xlabel('External CO_2 concentration (\muM)')
ylabel('Carboxylation Flux (units?)')
title('Carboxylation Flux In Different Models')
legend('Full CCM Carboxylation', 'Full CCM Oxygenation', ...
    'C3 RuBisCO Carboxylation', 'C3 RuBisCO Oxygenation', ...
    'Specific RuBisCO Carboxylation', 'Specific RuBisCO Oxygenation');
legend('boxoff');

% Fraction oxygenation as a function of CO2. 
frac_oxygenation_ccm = OratewC./CratewO;
frac_oxygenation_just_c3_rbc = OratewC_just_c3_rbc./CratewO_just_c3_rbc;
frac_oxygenation_just_specific_rbc = OratewC_just_specific_rbc./CratewO_just_specific_rbc;

figure(222);
loglog(CO2extv, frac_oxygenation_ccm, 'b');
hold on;
plot(CO2extv, frac_oxygenation_just_c3_rbc, 'g');
plot(CO2extv, frac_oxygenation_just_specific_rbc, 'r');
title('Ratio of Oxygenations/Carboxylations');
xlabel('external pH');
ylabel('Oxygenations/Carboxylations')
legend('CCM', 'C3 RuBisCO alone', 'Specific RuBisCO alone');

% Total cost a function of CO2.
figure(111);

plot(CO2extv, total_cost_just_c3_rbc_h, 'g');
hold on
plot(CO2extv, total_cost_just_specific_rbc_h, 'b');
fill([CO2extv fliplr(CO2extv)], ...
     [max_total_ccm_cost fliplr(min_total_ccm_cost)], 'm', ...
     'FaceAlpha', 0.3);
plot([15 15], [10 500], '--m');
axis([1e-4 1000 10 500]);
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
xlabel('External CO_2 concentration (\muM)');
ylabel('Energetic Cost H^+ / (CO_2 fixed)');
title('Total cost of fixation with CCM vs selective RuBisCO');
legend('C3 RuBisCO alone', 'Specific RuBisCO alone', 'Full CCM Range', ...
    'Atmospheric CO_2');
legend('boxoff');

