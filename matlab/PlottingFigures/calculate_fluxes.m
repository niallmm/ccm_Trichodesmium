function [fluxes] = calculate_fluxes(p, sol)
% calculate the fluxes from either analytic or numerical solution
% need p to be the parameterfile for the simultaion
% need sol to be the concentration solution for hte simulation
% sol.h_csome_uM, sol.h_cyto_uM, sol.c_csome_uM, sol.c_cyto_uM
% are the concetration values required to calculate all the fluxes
C = sol.c_csome_uM;
% carboxylation and oxygenation
fluxes.VO = p.VmaxCsome*p.KO/(p.Km*p.S_sat);
fluxes.CratewO_pm = p.Vmax*C./(C+p.Km*(1+p.O/p.KO))*p.Vcsome*1e3;
fluxes.CratewO_um = p.Vmax*C./(C+p.Km*(1+p.O/p.KO))*p.Vcsome*1e-3; % convert from uM*cm^3 to umoles
fluxes.OratewC_pm = fluxes.VO*p.O./(p.O+p.KO*(1+C/p.Km))*p.Vcsome*1e3;
fluxes.OratewC_um = fluxes.VO*p.O./(p.O+p.KO*(1+C/p.Km))*p.Vcsome*1e-3;
% error rate of oxygenation to total RuBisCO activity.
fluxes.error = fluxes.OratewC_pm/(fluxes.CratewO_pm + fluxes.OratewC_pm);
% HCO3- active flux into the cell
fluxes.Hin_pm = p.jc*p.Hout*p.SAcell*1e3;
fluxes.Hin_um = p.jc*p.Hout*p.SAcell*1e-3;

if isfield(sol, 'c_cyto_uM')
    % CO2 leakage across carboxysome
    fluxes.Ccsomeleak_um = -p.kcC*(sol.c_csome_uM - sol.c_cyto_rad_uM(1))*(4*pi*p.Rc^2)*1e-3;
    fluxes.Ccsomeleak_pm = fluxes.Ccsomeleak_um*1e6;
    % CO2 conversion
    fluxes.Cconv_pm = p.alpha*sol.h_cyto_uM*p.SAcell*1e3;
    fluxes.Cconv_um = p.alpha*sol.h_cyto_uM*p.SAcell*1e-3;
    % Net HCO3- leackage out of cell
    fluxes.Hleak_pm = (p.kmH_out*p.Hout - p.kmH_in*sol.h_cyto_uM)*p.SAcell*1e3;
    fluxes.Hleak_um = (p.kmH_out*p.Hout - p.kmH_in*sol.h_cyto_uM)*p.SAcell*1e-3;
    % net CO2 leakage out of cell
    fluxes.Cleak_um = p.kmC*(p.Cout - sol.c_cyto_uM)*p.SAcell*1e-3;
    fluxes.Cleak_pm = p.kmC*(p.Cout - sol.c_cyto_uM)*p.SAcell*1e3;
    
    % production of OH due to HCO3- conversion
    fluxes.OHrate_um = -p.kcH*(sol.h_csome_uM - sol.h_cyto_uM)*(4*pi*p.Rc^2)*1e-3;
    fluxes.OHrate_pm = -p.kcH*(sol.h_csome_uM - sol.h_cyto_uM)*(4*pi*p.Rc^2)*1e3;
    fluxes.OHrate_alt_um = p.kcC*(sol.c_csome_uM - sol.c_cyto_uM)*(4*pi*p.Rc^2)*1e-3...
        + fluxes.CratewO_um;
    fluxes.OHrate_alt_pm = p.kcC*(sol.c_csome_uM - sol.c_cyto_uM)*(4*pi*p.Rc^2)*1e3 ...
        + fluxes.CratewO_pm;
    
    
else  % the csome is actually the whole cell.
    % CO2 conversion
    fluxes.Cconv_pm = p.alpha*sol.h_csome_uM*p.SAcell*1e3;
    fluxes.Cconv_um = p.alpha*sol.h_csome_uM*p.SAcell*1e-3;
    % Net HCO3- leackage out of cell
    fluxes.Hleak_pm = (p.kmH_out*p.Hout - p.kmH_in*sol.h_csome_uM)*p.SAcell*1e3;
    fluxes.Hleak_um = (p.kmH_out*p.Hout - p.kmH_in*sol.h_csome_uM)*p.SAcell*1e-3;
    % net CO2 leakage out of cell
    fluxes.Cleak_um = p.kmC*(p.Cout - sol.c_csome_uM)*p.SAcell*1e-3;
    fluxes.Cleak_pm = p.kmC*(p.Cout - sol.c_csome_uM)*p.SAcell*1e3;
    
end
end
