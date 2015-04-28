
% Bicarbonate flux vs facilitated uptake  vs scavenging

phoutv = linspace(6, 8.5, 20);

addpath('/Users/niallmangan/GitHub/ccm/matlab')
p = CCMParams_Csome;
p.pH = 8;
p.jc = 2e-4;
p.alpha = p.jc;
p.kRub = 11.6; % rxns/s maximum reaction rate at single active site
p.Km_8 = 340;    % half max reaction rate of RuBisCO, uM
p.S_sat = 43;  % specificity ratio
p.KO = 972;    % uM
p.Ci_tot = 15; 

exec = FullCCMModelExecutor(p);
res = exec.RunAnalytical();

fignum = 11;
 klist =[1e-4 1e-2 1];
% for i = 1:100:length(HoutX)
for m = 1:length(klist)
    p.k = klist(m);
    for i = 1:length(phoutv)
          p.pH_out = phoutv(i);
          exec = FullCCMModelExecutor(p);
          res = exec.RunAnalytical();
        % all converted to picomoles
        CoutX(i) = p.Cout;
        HoutX(i) = p.Hout;
        Htransport(i) = res.Hin_pm;
        Cfacilitateduptake(i) = p.alpha*p.kmC*p.G*CoutX(i)*1e6*4*pi*p.Rb^2/...
            ((p.alpha+p.kmC)*p.G+p.D/p.Rb^2);
        Cscavenging(i) = p.alpha*res.c_csome_uM*1e6*4*pi*p.Rb^2*(1-(p.alpha+p.kmC)*p.G/((p.alpha+p.kmC)*p.G + p.D/p.Rb^2));
        
        
        
    end
    

   
    if m == 1
        k_plotting = p.k
        figure(fignum)
        semilogy(phoutv, Htransport./HoutX, 'b')
        hold on
        plot(phoutv, Cfacilitateduptake./CoutX, 'r')
    end

    figure(fignum)
    semilogy(phoutv, Cscavenging./p.Ci_tot, '-.r')

    
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
xlabel('external pH')
ylabel('Flux of [C_i] at cell membrane [picomoles/s]')
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
haxes2 = axes('Position',haxes1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Yscale', 'log', ...
    'Color','none');
hold on
axis([CoutX(end)/(CoutX(end)+HoutX(end)) CoutX(1)/(CoutX(1)+HoutX(1)) 1e-9 1e-2])
set(haxes2,'ytick',[])
set(haxes2,'yticklabel',[])
set(haxes2, 'XDir','Reverse')
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