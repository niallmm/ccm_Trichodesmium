function total_cost_h = totalProtonCost_CCM(Hin, CratewO, ...
    OratewC, OHrateCA, BicarbonateUptakeCost)    
    cost_fix_and_recover = fixAndRecover_ProtonCost(CratewO, OratewC);
    
    % Note: Luke makes an excellent point that there is a Maxwell's demon
    % here. We assume transport costs 1 H+ per HCO3- regardless of the
    % concentration of HCO3- inside and outside the cell. If H+ has a
    % roughly constant potential due to the pH difference across the cell
    % membrane, then this doesn't make sense. The cost of transport should
    % go up as the extracellular concentration of CO2 goes down. 
    cost_h = BicarbonateUptakeCost * Hin./CratewO;
    cost_pH_maintenance = (Hin - OHrateCA)./CratewO;
    
    total_cost_h = cost_fix_and_recover + cost_h + cost_pH_maintenance;
end

