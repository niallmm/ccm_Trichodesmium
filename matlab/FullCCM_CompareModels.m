%==========================================================================
% Compares the results of the analytical and numerical models in the Full
% CCM case
%==========================================================================
clear all
%  close all
 kcvec =logspace(-12,1,20);
%  jcvec = logspace(-6,1, 20);
%  jcvec= 1e-4;
 
for ii = 1:length(kcvec)
% for ii = 1:length(jcvec)
% parameters for calculation.
p = CCMParams_Csome();
p.pHoff = 1;
% p.jc = jcvec(ii);
p.kcC = kcvec(ii);
p.kcH = kcvec(ii);
% p.kcC = 1;
% p.kcH = 1;
p.pH = 8;
p.pH_out = 8;

p.jc = 1e-6;
%  p.NCA = 60;
% p.NCA = 200;
 p.Cout = 15;
 p.I_out = p.I_in;
%    p.Cout = 0.24;


p.alpha = 1e-5;
%  p.Rc = 1e-6;


% % Optimal jc for the no csome case using analytic solution
% Hmax = 30000; % 30 mM in uM
% jc_opt= p.CalcOptimalJc(Hmax); 
% p.jc = jc_opt;
% 
% % Numerical solution 
executor = FullCCMModelExecutor(p);
num_res = executor.RunNumerical();
h_csome_mMn(ii) = num_res.h_csome_mM;
c_csome_mMn(ii) = num_res.c_csome_mM;
r_cell = num_res.r_cell;
h_cyto_mMn = num_res.h_cyto_rad_mM;
c_cyto_mMn = num_res.c_cyto_rad_mM;

% Analytical solution from the no csome case
ana_res = executor.RunAnalytical();
r_cella = ana_res.r;
h_cyto_mMa = ana_res.h_cyto_rad_uM*1e-3;
c_cyto_mMa = ana_res.c_cyto_rad_uM*1e-3;
h_csome_mMa(ii) = ana_res.h_csome_mM;
c_csome_mMa(ii) = ana_res.c_csome_mM;
c_csome_sat(ii) = ana_res.CCAsat0*1e-3;


% non dimenstional analytic solution
c_nondim(ii) = p.Kca*(p.kappa*p.epsilon_h+p.epsilon_c+p.Vmax*p.Rc^2/p.Kca/3/p.D)./(p.beta_h*p.kappa*p.gamma + p.beta_c)*1e-3;

end
% figure(1)
% loglog(jcvec, c_csome_mMa, 'r')
% hold on
% loglog(jcvec, h_csome_mMa, 'b')
% loglog(jcvec, c_csome_mMn, 'or')
% loglog(jcvec, h_csome_mMn, 'ob')
% loglog(jcvec, c_nondim, '--k')
% xlabel('HCO_3^- transport')
% ylabel('CO_2 and HCO_3^- concentration in mM')

figure(1)
loglog(kcvec, c_csome_mMa, 'r')
hold on
loglog(kcvec, h_csome_mMa, 'b')
loglog(kcvec, c_csome_mMn, 'or')
loglog(kcvec, h_csome_mMn, 'ob')
xlabel('carboxysome permeability')
ylabel('CO_2 and HCO_3^- concentration in mM')
% 
% figure(5)
% semilogy(r_cell, c_cyto_mMn, 'or')
% hold on
% plot(r_cella, c_cyto_mMa , '-r')
% plot(r_cell, h_cyto_mMn, 'ob')
% plot(r_cella, h_cyto_mMa, '-b')
%  line([0 p.Rc],[c_csome_mMn(end) c_csome_mMn(end)], 'Color','r', 'Marker', 'o') 
%  line([0. p.Rc], [c_csome_mMa(end) c_csome_mMa(end)], 'Color', 'r', 'Marker', '*')
%  line([0 p.Rc], [h_csome_mMn(end) h_csome_mMn(end)], 'Color', 'b', 'Marker', 'o')
%  line([0. p.Rc], [h_csome_mMa(end) h_csome_mMa(end)], 'Color', 'b', 'Marker', '*')
% xlabel('Cell radius (cm)')
% ylabel('Concentration (mM)')
%  line([p.Rc p.Rc], [1e-6 1e6], 'Color', [0.5 0.5 0.5], 'LineStyle', '-.', 'LineWidth', 3)
%  xlim([0 5e-5])
%  ylim([1e-5 100])
% set(gca,'XTick',[0 1e-5 2e-5 3e-5 4e-5 5e-5])
% % set(gca,'YTick',[1e-5  1e-3  1e-1  10])
% box off