function total_cost_h = fixAndRecover_ProtonCost(CratewO, OratewC)
    % This is the calculated cost of fixing CO2, recovering oxygenation
    % products and re-fixing them. The calculation itself is described in
    % detail in the SI of the paper. Note that this calculation accounts
    % for the inherently recursive nature of this process: photorespiration
    % decarboxylates, producing CO2 which necessitates another fixation,
    % which may lead to another round of oxygenation etc.
    
    % Note that we are assuming here that the rate of decarboxylation
    % through photorespitation is not sufficient to change the cellular
    % concentration of CO2. Since we are here handling the case where there
    % is no CCM, the ratio of oxygenation to carboxylation rates may be as
    % high as 30% (Sharkey, 1988). If the absolute rate of carboxylation is
    % on the order of 1e-20 mol/cell/s (as produced by the model here and
    % measured in Hopkinson et al, 2014) then the absolute rate of
    % photorespiration is at most 3.3e-21 mol/cell/s, producing CO2 at the
    % same rate. 
    
    % If that rate of CO2 production was trapped in a model cell with
    % volume 4/3*pi*(5e-5 cm)^3 = 5.2e-16 L = 0.5 fL, it would produce an
    % additional 6 uM CO2 / cell / s. Since the membrane permeability to
    % CO2 is here taken to be 0.3 cm/s that accumulated 6 uM of cytosolic
    % CO2 will cross the membrane diffusively at a rate 0.3 cm/s * 6 uM =
    % 1.8e-9 mol/cm^2/s. Multiplying by the cell surface area, which is
    % 4*pi*(5e-5cm)^2 = 3.14e-8 cm^2, we get 1.8e-9*3.14e-8 = 5.7e-17
    % mol/cell/s. So even in the limiting case of high photorespiration, 
    % the rate at which CO2 leaks through the membrane is 4 orders of
    % magnitude faster than the rate at which it is produced through
    % photorespitation. As such, it is reasonable to neglect the effects 
    % of photorespiratory CO2 production on the cytosolic CO2 concentration.
    ratio = OratewC ./ CratewO;
    total_cost_h = (33.0 + (47.0 * ratio)) ./ (1.0 - (0.5 * ratio));
end

