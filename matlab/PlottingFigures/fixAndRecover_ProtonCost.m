function total_cost_h = fixAndRecover_ProtonCost(CratewO, OratewC)
    % This is the calculated cost of fixing CO2, recovering oxygenation
    % products and re-fixing them. The calculation itself is described in
    % detail in the SI of the paper. Note that this calculation accounts
    % for the inherently recursive nature of this process: photorespiration
    % decarboxylates, producing CO2 which necessitates another fixation,
    % which may lead to another round of oxygenation etc.
    ratio = OratewC ./ CratewO;
%     total_cost_h = (32.67 + (49.0 * ratio)) ./ (1.0 - (0.5 * ratio));
    total_cost_h = (36.67 + (56*ratio))./ (1.0 - (0.5*ratio)); % new number for H+/NADPH from Raven
end

