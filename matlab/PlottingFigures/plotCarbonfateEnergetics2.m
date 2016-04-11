%Run calcCarbonfate_vs_pH.m first

all_ones = ones(size(pH));

% cost of CBB per carboxylation
% cost_fixation = 32.67;
cost_fixation = 36.67; % with Raven update
% cost of photorespiration + CO2 recovery per oxygenation
% cost_per_2pg = 65.34;
cost_per_2pg = 74.34; % with Raven update
% Cost of transport in the 3 cases (per fixation)
cost_h = Hin./CratewO;



%Cost of maintaining pH balance with a CCM 
% - only account for pH maintanence needed from CCM activity
% -- 1 ATP per HCO3- transported (need to export 1 H+?)
% -- gain 1 OH group (or loose one H+) from each HCO3- converted to CO2
cost_pH_maintenance = (Hin - OHrateCA)./CratewO;


% use the oxygenations and carboxylation rates calculated in the analytic
% and model execut`or files (calc_Carbonfate_vs_pH.m calls them)
frac_oxygenation_ccm = OratewC./CratewO;

% Cost of 2PG recovery on a per-fixation basis in each case.
% cost_2pg_ccm_h2 = cost_per_2pg * frac_oxygenation_ccm;
cost_2pg_ccm_h = fixAndRecover_ProtonCost(CratewO, OratewC)-cost_fixation;

% Total cost of the system in H+/fixation for each case.
% costpertransport = 4;
total_cost_ccm_h = totalProtonCost_CCM(Hin, CratewO, OratewC, OHrateCA, costpertransport);


%%
% =========================================================================
% Plot cost of transport, CCB cycle and photorespiration per carbon fixation
% =========================================================================
figure(1)
% for model with full CCM


% plot cost of transport/CO2 fixation with pH
semilogy(pH, cost_h,'Color','r', 'LineWidth', 3);
axis([pH(1) pH(end) 1e-1 2e4])
hold on
ylabel('Energetic Cost H^+ / (CO_2 fixed)')
xlabel('Cytoplasmic pH')
text(6, 5e2, 'HCO3^- transport cost', 'Color', 'r')
% plot cost of fixation/CO2 (flat)
plot([pH(1) pH(end)], [cost_fixation cost_fixation], '--b');
text(3, 20, 'Calvin Cycle Fixation', 'Color', 'b');
% plot cost of photorespiration
semilogy(pH, cost_2pg_ccm_h, 'Color','y', 'LineWidth', 3);
text(6, 1, 'Cost of photorespiration', 'Color', 'k')
%plot cost of pH maintenance and transport
semilogy(pH, cost_pH_maintenance, '--c', 'LineWidth',3)
% total cost
semilogy(pH, total_cost_ccm_h, 'k', 'Linewidth', 3)
title('Full CCM')
% 
line([8.39 8.39],[0.01 1e3])
line([7.33 7.33],[0.01 1e3])
line([8.28 8.28],[0.01 1e3])
line([8.49 8.49],[0.01, 1e3])
line([7.17 7.17], [0.01, 1e3])
line([7.49 7.49],[0.01, 1e3])
drawnow

