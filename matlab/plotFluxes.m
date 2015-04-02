% =========================================================================
% Plot CO2 and HCO3^- fluxes as a function of pH
% =========================================================================

    figure
    plot(pH, abs(Hin), 'r')
    hold on
    plot(pH, abs(CratewO), 'b')
    plot(pH, abs(OratewC), 'c')
    plot(pH, abs(Hleak), 'k')
    plot(pH, abs(Cleak), 'g')
    xlabel('pH')
    ylabel('inorganic carbon fluxes (picomoles/s)')
    legend('HCO_3^- transport', 'carboxylation', 'oxygenation','HCO_3^- leakage', 'CO_2 leakage')

