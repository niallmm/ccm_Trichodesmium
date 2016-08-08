% =========================================================================
% Plot CO2 and HCO3^- fluxes as a function of pH
% =========================================================================
% run calcCarbonfate_vs_pH first. 

% rate consants with  S = 10, approximating I = 0.2M.
    kp1 = 3.711e-2; %1/s
    km1 = 3.649e4; %M/s
    kp4 = 1780; %M/s
    km4 = 5.26e-5; %1/s
    
    kd = km1*10.^-pH + km4;
    kh = kp1 + kp4*10.^(pH-14);
    
    Hloss = (kd.*Hmax - kh.*Ccyto)*ccm_params.Vcell*1e3; % uM/s -> umol/s
    
    figure
    semilogy(pH, abs(Hin), 'r')
    hold on
    
    plot(pH, abs(Hleak), 'k')
    plot(pH, abs(Cleak), 'g')
    semilogy(pH, abs(CratewO), 'b')
    plot(pH, abs(Hloss), 'c')
    xlabel('cytosolic pH')
    ylabel('inorganic carbon fluxes (picomoles/s)')
    legend('HCO_3^- transport', 'HCO_3^- leakage', 'CO_2 leakage', 'carboxylation', 'HCO_3^- to CO_2 natural conversion')
    legend('boxoff')
    
    figure
    plot(pH, Hloss./Hin, 'k')
    xlabel('cytosolic pH')
    ylabel('ratio of cytosolic dehydration \newline rate to HCO_3^- transport')
    
	figure
    plot(pH, abs(Hloss./Cleak))
    xlabel('cytosolic pH')
    ylabel('ratio of cytosolic dehydration \newline rate to CO_2 leakage')

