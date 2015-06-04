
% Bicarbonate flux vs facilitated uptake  vs scavenging

phoutv = linspace(6, 8.5, 20);

addpath('/Users/niallmangan/GitHub/ccm/matlab')
p = CCMParams_Csome;
p.pH = 8;
p.jc = 2e-3;
p.alpha = p.jc;
p.kRub = 11.6; % rxns/s maximum reaction rate at single active site
p.Km_8 = 340;    % half max reaction rate of RuBisCO, uM
p.S_sat = 43;  % specificity ratio
p.KO = 972;    % uM
p.Ci_tot = 15000;

exec = FullCCMModelExecutor(p);
res = exec.RunAnalytical();

fignum = 10;
 klist =[1e-4 1e-2 1];
% for i = 1:100:length(HoutX)
for m = 1:length(klist)
    p.k = klist(m);
    for i = 1:length(phoutv)
          p.pH_out = phoutv(i);
          exec = FullCCMModelExecutor(p);
          res = exec.RunAnalytical();
        % all converted to picomoles
        Htransport(i) = res.Hin_pm;
        Cfacilitateduptake(i) = p.alpha*p.kmC*p.G*CoutX(i)*1e6*4*pi*p.Rb^2/...
            ((p.alpha+p.kmC)*p.G+p.D/p.Rb^2);
        Cscavenging(i) = p.alpha*res.c_csome_uM*1e6*4*pi*p.Rb^2*(1-(p.alpha+p.kmC)*p.G/((p.alpha+p.kmC)*p.G + p.D/p.Rb^2));
        
        
        
    end
    

   
    if m == 1
        k_plotting = p.k
        figure(fignum)
        semilogy(CoutX./(HoutX +CoutX), Htransport, 'b')
        hold on
        plot(CoutX./(HoutX +CoutX), Cfacilitateduptake, 'r')
    end

    figure(fignum)
    semilogy(CoutX./(HoutX +CoutX), Cscavenging, '-.r')

    
    if m==1
        k_legend = klist(1)
        figure(fignum)
        hleg = legend('HCO_3^- Active Transport', 'Facilitated CO_2 uptake','CO_2 Scavenging', 'Location', 'Best');
        legend('boxoff')
    end
    
    
end

figure(fignum)

% plot(CoutX./(HoutX + CoutX), Cleakage, 'k')
% plot(CoutX./(HoutX + CoutX), Hleakage, '--k')
% plot(CoutX./(HoutX + CoutX), Cleakage+Hleakage, 'g')
xlabel('Proportion of external C_i composed of [CO_2]: \newline [CO_2]/([HCO_3^-] + [CO_2])')
ylabel('Flux of [HCO_3^-] at cell membrane [picomoles/s]')

% figure(11)
% plot(CoutX./(HoutX +CoutX), HcytoRbstore, 'b')
% hold on
% plot(CoutX./(HoutX +CoutX), CcytoRbstore, 'r')
% xlabel('Proportion of external C_i composed of [CO_2]')
% ylabel('Concentration at cell membrane [mM]')
% figure(5)
% hold on
% plot(CoutX./(HoutX +CoutX), Htransport+Cfacilitateduptake+Cscavenging, 'b')
%



% rc = linspace(0,Rc, 1e2);
% rb = linspace(Rc, Rb, 1e3);
%
%
% Hcsome = Ccsomep*Vca*Kba/(Vba*Kca);
%
% Hcyto = ((jc+kmH)*HoutX + alpha*CcytoRb - kmH*Hcsome)*...
%     (D/(k*Rc^2)+1/Rc -1./rb)/(kmH*G + D/Rb^2) + Hcsome;