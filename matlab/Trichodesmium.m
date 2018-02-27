% Simulations for "Inorganic C and pH dependence of photosynthetic rates of
% Trichodesmium." experiments and paper by Tobias G. Boatman (1), Tracy Lawson (1), 
% and Richard J. Geider (1). Simulations run by Niall Mangan (2)
% parameters updated for Trichodesmium
% (1) School of Biological Sciences, University of Essex, Wivenhoe Park, Colchester, Essex, CO4 3SQ, UK
% (2) Department of Engineering Sciences and Applied Mathematics, Northwestern University, Evanston, IL, USA 60208.
clear all
close all
changeplot
addpath('PlottingFigures')
addpath(fileparts(pwd))

% Params for cell with full CCM.
ccm_params = CCMParams_Csome;
nocsome_params = CCMParams_NoCsome;

% update RuBisCO params
NRub = 2000;
baseV = ccm_params.Vcsome;
NCA = 100;
baseSA = 4*pi*ccm_params.Rc^2;
ccm_params.kRub = 1.92;      % rxns/s maximum reaction rate at single active site
ccm_params.Km_meas = 145;    % half max concentration for CO2, uM
ccm_params.S_sat = 45;      % specificity ratio
ccm_params.KO_meas = 600;     % half max concentration for O2, uM

% update geometry
ccm_params.Rb = 3.058e-4; %cm
ccm_params.Rc = 1.5e-5; %cm

% scale up number of RuBisCO using csome volume:
ccm_params.NRub= NRub*ccm_params.Vcsome/baseV;
% scale up number of CA using csome surface area:
ccm_params.NCA = NCA*(4*pi*ccm_params.Rc^2)/baseSA;

% set internal pH
ccm_params.pH = 8.3; % internal pH
ccm_params.pH_csome = 8.3; % internal pH
ccm_params.pHoff = 0;
% salt water system:
ccm_params.salt = 1;

% pH of 7.5 to 8.5, and external CO2 range of about 0.003 to 0.03 mmol/L-1.
% set external pH
ccm_params.pH_out = 7.5;
ccm_params.Cout = 9; % 3 uM to 30 uM
ccm_params.Hout = 1900;%ccm_params.Cout*10^(-ccm_params.pKa_eff_out+ccm_params.pH_out);

% hack the equilibrium constant in the numerical simulation
% this makes Keq the correct value for the internal pH
Keq = ccm_params.Keq;
ccm_params.Kca = ccm_params.Vca*ccm_params.Kba/(ccm_params.Vba*ccm_params.Keq);

% set carboxysome permeability to that in Mangan et al 2016
ccm_params.kcC = 3e-5;
ccm_params.kcH = 3e-5;



% set jc for desired internal HCO3- concentration
jc_opt = 3e-7; % for Hmax=  30 mM;
ccm_params.jc = jc_opt*0.8;
ccm_params.alpha = 0.2*jc_opt; %attribute some level activity to 
% CO2->HCO3- conversion to sustain internal HCO3- levels

% run problem
exec = FullCCMModelExecutor(ccm_params);
num = exec.RunNumerical();
% calculate fluxes from model output
[fluxes] = calculate_fluxes(ccm_params, num);

%plot concentration gradients in cell for these baseline values.
figure(1)
plot(num.r*ccm_params.Rc, num.c_mM(end,:))
hold on
plot(num.r_cell, num.c_cyto_rad_mM)
plot(num.r*ccm_params.Rc, num.h_mM(end,:))
plot(num.r_cell, num.h_cyto_rad_mM)
%% fixed pH, varying totIC, %CO2 from 3 uM to 30 uM
ccm_params.pH_out = 8.15;
CO2extvary = linspace(0.05, 30,30);

for ii = 1:length(CO2extvary)
    ccm_params.Cout = CO2extvary(ii);
    ccm_params.Hout = ccm_params.Cout*10^(-ccm_params.pKa_eff_out+ccm_params.pH_out);
    Hout1(ii) = ccm_params.Hout;
    exec = FullCCMModelExecutor(ccm_params);
    num = exec.RunNumerical();
    num1(ii) = num;
    C_csome1(ii) = num.c_csome_mM;
    H_csome1(ii) = num.h_csome_mM;
    fluxes(ii) = calculate_fluxes(ccm_params, num);
end
% 
figure
semilogy(CO2extvary, abs([fluxes.Hin_um]), 'Color',newcolor(8,8))
hold on
plot(CO2extvary, abs([fluxes.Hleak_um]), 'Color',newcolor(4,8))
plot(CO2extvary, abs([fluxes.Cleak_um]), 'Color',newcolor(3,8))
plot(CO2extvary, abs([fluxes.CratewO_um]),'Color',newcolor(1,8))
xlabel('external CO_2 [mM]')
ylabel('fluxes [uM/s] for model cell')
 legend('HCO_3^- transport', 'HCO_3^- leakage', 'CO_2 leakage', 'carboxylation')
axis([min(CO2extvary) max(CO2extvary) 1e-15 1e-11])
 title('varying CO_2 fixed pH')
%% pH vary, CO2 vary
pHextvary = linspace(7.52, 8.54, 20);
ccm_params.Hout = 1900;
for ii = 1:length(pHextvary)
    ccm_params.pH_out = pHextvary(ii);
    CO2extvary1(ii) = ccm_params.Hout*10^(ccm_params.pKa_eff_out-ccm_params.pH_out);
    ccm_params.Cout = CO2extvary1(ii);
    HCO3out(ii) = ccm_params.Hout;
    exec = FullCCMModelExecutor(ccm_params);
    num = exec.RunNumerical();
    num2(ii) = num;
    C_csome_pH_vary(ii) = num.c_csome_mM;
    H_csome_pH_vary(ii) = num.h_csome_mM;
    fluxes_ph_vary(ii) = calculate_fluxes(ccm_params, num);
end
%
figure
semilogy(CO2extvary1, abs([fluxes_ph_vary.Hin_um]), 'Color',newcolor(8,8))
hold on
plot(CO2extvary1, abs([fluxes_ph_vary.Hleak_um]), 'Color',newcolor(4,8))
plot(CO2extvary1, abs([fluxes_ph_vary.Cleak_um]), 'Color',newcolor(3,8))
plot(CO2extvary1, abs([fluxes_ph_vary.CratewO_um]),'Color',newcolor(1,8))
xlabel('external CO_2 [mM]')
ylabel('fluxes [uM/s] for model cell')
 legend('HCO_3^- transport', 'HCO_3^- leakage', 'CO_2 leakage', 'carboxylation')
axis([min(CO2extvary1) max(CO2extvary1) 1e-14 5e-11])
 title('varying CO_2 varying pH')
 
 %% CO2 fixation for both experiments
 figure
plot(CO2extvary1, abs([fluxes_ph_vary.CratewO_um]),'o','Color',newcolor(2,8))
hold on
plot(CO2extvary, abs([fluxes.CratewO_um]), 'o','Color',newcolor(1,8),'MarkerFaceColor', newcolor(1,8)) 
 xlabel('external CO_2 [uM]')
 ylabel('carboxylation rate [uM/s]')
 legend('pH varies', 'HCO_3^- varies')
 legend boxoff
 %% Plot as a function of HCO3-
 figure
 plot(Hout1, abs([fluxes.CratewO_um]), 'o','Color',newcolor(1,8),'MarkerFaceColor', newcolor(1,8)) 
 hold on
 plot(HCO3out,  abs([fluxes_ph_vary.CratewO_um]),'o','Color',newcolor(2,8))
 xlabel('external HCO_3^- [uM]')
 ylabel('carboxylation rate [uM/s]')
 legend('HCO_3^- varies', 'pH varies')
 legend boxoff
 %% plot as a function of HCO3- uptake - H leakage
 figure
 plot([fluxes.Hin_um]+[fluxes.Hleak_um], abs([fluxes.CratewO_um]), 'o','Color',newcolor(1,8),'MarkerFaceColor', newcolor(1,8)) 
 hold on
 plot([fluxes_ph_vary.Hin_um]+[fluxes_ph_vary.Hleak_um],  abs([fluxes_ph_vary.CratewO_um]),'o','Color',newcolor(2,8))
 xlabel('HCO_3^- uptake - HCO_3^- leakage [uM/s]')
 ylabel('carboxylation rate [uM/s]')
  legend('HCO_3^- varies', 'pH varies')
  legend boxoff
 
 
%% compare numerical and analytic solutions for varying uptake rates
jcvec= logspace(-8, -3, 20);

for ii = 1:length(jcvec)
    ccm_params.jc = jcvec(ii);
    exec = FullCCMModelExecutor(ccm_params);
     res = exec.RunAnalytical();
    num = exec.RunNumerical();
    
   C_csome(ii) = num.c_csome_mM;
    H_csome(ii) = num.h_csome_mM;
     
     C_csome_ana(ii) = res.c_csome_mM;
     H_csome_ana(ii) = res.h_csome_mM;
    
end

figure
loglog(jcvec, C_csome,'o', 'MarkerFaceColor',newcolor(1,2))
hold on
loglog(jcvec, H_csome,'o', 'MarkerFaceColor',newcolor(2,2))
loglog(jcvec, C_csome_ana,'-', 'MarkerFaceColor',newcolor(1,2))
loglog(jcvec, H_csome_ana,'-', 'MarkerFaceColor',newcolor(2,2))

xlabel('HCO_3^- transport [cm/s]')
ylabel('CO_2 concentration [mM]')
%% compare analytic and numerical solutions for varying carboxysome permeabilities
kcvec= logspace(-7, 2, 30);
ccm_params.jc = 1e-4;
for ii = 1:length(kcvec)
    ccm_params.kcC = kcvec(ii);
    ccm_params.kcH = kcvec(ii);
    exec = FullCCMModelExecutor(ccm_params);
     res = exec.RunAnalytical();
    num = exec.RunNumerical();
    
    C_csome(ii) = num.c_csome_mM;
    H_csome(ii) = num.h_csome_mM;
    
    C_csome_ana(ii) = res.c_csome_mM;
    H_csome_ana(ii) = res.h_csome_mM;
    
end

figure
loglog(kcvec, C_csome,'o', 'MarkerFaceColor',newcolor(1,2))
hold on
loglog(kcvec, H_csome,'o', 'MarkerFaceColor',newcolor(2,2))
loglog(kcvec, C_csome_ana, '-', 'MarkerFaceColor', newcolor(1,2))
loglog(kcvec, H_csome_ana, '-', 'MarkerFaceColor', newcolor(2,2))
xlabel('carboxysome permeability [cm/s]')
ylabel('CO_2 concentration [mM]')