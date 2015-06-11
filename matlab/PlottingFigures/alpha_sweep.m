
% Bicarbonate flux vs facilitated uptake  vs scavenging

phoutv = linspace(6, 8.5, 20);

addpath('/Users/niallmangan/GitHub/ccm/matlab')
p = CCMParams_Csome;
p.pH = 8;
%p.jc = 2e-4;
p.jc = 0;
alphasweep = logspace(-8, -2,20);

% p.alpha= 0;
p.kRub = 11.6; % rxns/s maximum reaction rate at single active site
p.Km_8 = 340;    % half max reaction rate of RuBisCO, uM
p.S_sat = 43;  % specificity ratio
p.KO = 972;    % uM
% p.Ci_tot = 15; 

exec = FullCCMModelExecutor(p);
res = exec.RunAnalytical();

fignum = 12;
for i = 1:length(alphasweep)
p.alpha = alphasweep(i);
p.jc = p.alpha;
          exec = FullCCMModelExecutor(p);
          res = exec.RunAnalytical();
        % all converted to picomoles
        CoutX(i) = p.Cout;
        HoutX(i) = p.Hout;
        Hcyto(i) = res.h_cyto_mM;
        Ccyto(i) = res.c_cyto_mM;
        Hcsome(i) = res.h_csome_mM;
        Ccsome(i) = res.c_csome_mM;
        Cleak(i) = res.Cleak_pm;
        Hleak(i) = res.Hleak_pm;
        Cfix(i) = res.CratewO_pm;
        Htransport(i) = res.Hin_pm;
        Cfacilitateduptake(i) = p.alpha*p.kmC*p.GC*(CoutX(i)-res.c_csome_uM)*1e3*4*pi*p.Rb^2/...
            ((p.alpha+p.kmC)*p.GC+p.D/p.Rb^2);
        Cscavenging(i) = p.alpha*res.c_csome_uM*1e3*4*pi*p.Rb^2*(1-(p.alpha+p.kmC)*p.GC/((p.alpha+p.kmC)*p.GC + p.D/p.Rb^2));
        
        Cconversion(i) = p.alpha*res.c_cyto_uM*1e3*4*pi*p.Rb^2;
        
        



%    
%    % if m == 1
%         k_plotting = p.kcH
%         figure(fignum)
%         semilogy(alphasweep, Htransport, 'b')
%         hold on
%          plot(alphasweep, Cconversion, 'k')
% %         figure(223)
%        plot(alphasweep, abs(Cfacilitateduptake), 'r')
%         hold on
%         
%  %   end
% 
% %    figure(fignum)
%     semilogy(phoutv, Cscavenging, '-.r')
% 
%     
%     if m==1f
%         k_legend = klist(1)
%         figure(fignum)
%    %     hleg = legend('HCO_3^- Active Transport', 'Facilitated CO_2 uptake','CO_2 Scavenging', 'Location', 'Best');
%    %     legend('boxoff')
%     end
    
    
end

figure(fignum)
   figure(fignum)
   loglog(alphasweep, Ccsome, 'r')
   hold on
   plot(alphasweep, Hcsome, 'b')
   plot(alphasweep, Hcyto, '--b')
   plot(alphasweep, Ccyto, '--r')
% plot(CoutX./(HoutX + CoutX), Cleakage, 'k')
% plot(CoutX./(HoutX + CoutX), Hleakage, '--k')
% plot(CoutX./(HoutX + CoutX), Cleakage+Hleakage, 'g')
xlabel('alpha')
ylabel('cytosolic C_i concentration')

figure(fignum+1)
loglog(alphasweep, Cleak, 'r')
hold on
plot(alphasweep, abs(Hleak), 'b')
% haxes1 = gca; % handle to axes
% haxes1_pos = get(haxes1,'Position'); % store position of first axes
% haxes2 = axes('Position',haxes1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','left',...
%     'Yscale', 'log', ...
%     'Color','none');
% hold on
% axis([CoutX(end)/(CoutX(end)+HoutX(end)) CoutX(1)/(CoutX(1)+HoutX(1)) 1e-9 1e-2])
% set(haxes2,'ytick',[])
% set(haxes2,'yticklabel',[])
% set(haxes2, 'XDir','Reverse')
% % figure(11)
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