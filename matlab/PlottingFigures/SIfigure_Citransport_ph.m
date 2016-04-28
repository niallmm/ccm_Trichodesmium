
% Bicarbonate flux vs facilitated uptake  vs scavenging
addpath(fileparts(pwd))
phoutv = linspace(6, 8.5, 20);

addpath('/Users/niallmangan/GitHub/ccm/matlab')
p = CCMParams_Csome;
p.pH = 8;
p.pH_out = 7;
kopt= 3e-5;
p.kcH = kopt;
p.kcC = p.kcH;

Hmax = 30000;   % Maximum cytoplasmic bicarbonate conc. in uM
jc_opt = p.CalcOptimalJc(Hmax);

exec = FullCCMModelExecutor(p);
res = exec.RunAnalytical();

fignum = 16;

for jj = 1:2
    if jj == 2
        p.alpha = jc_opt;
        p.jc = 0;
    else
        p.jc = jc_opt;
        p.alpha = 0;
    end
    for i = 1:length(phoutv)
        p.pH_out = phoutv(i);
        p.salt = 0;
        exec = FullCCMModelExecutor(p);
        res = exec.RunAnalytical();
        % all converted to picomoles
        CoutX(i) = p.Cout;
        HoutX(i) = p.Hout;
        Hcyto(i) = res.h_cyto_uM;
        Ccyto(i) = res.c_cyto_uM;
        Htransport(i) = res.Hin_pm;
        Cfacilitateduptake(i) = p.alpha*p.kmC*p.GC*(CoutX(i)-res.c_csome_uM)*1e3*4*pi*p.Rb^2/...
            ((p.alpha+p.kmC)*p.GC+p.D/p.Rb^2);
        Cscavenging(i) = p.alpha*res.c_csome_uM*1e3*4*pi*p.Rb^2*(1-(p.alpha+p.kmC)*p.GC/((p.alpha+p.kmC)*p.GC + p.D/p.Rb^2));
        
        Cconversion(i) = p.alpha*res.c_cyto_uM*1e3*4*pi*p.Rb^2;
        
        Cleak(i) = res.Cleak_pm; % positive going into cell (cytosol)
        Cleak_csome(i) = -p.kcC*(res.c_cyto_uM - res.c_csome_uM)*1e3*4*pi*p.Rc^2;
        % positive going into cytosol
        Hleak(i) = res.Hleak_pm;
        Cfix(i) = res.CratewO_pm;
        
        
    end
    
    
    p.pH_out = 8;
    p.salt = 1;
    exec = FullCCMModelExecutor(p);
    res = exec.RunAnalytical();
    % all converted to picomoles
    if jj == 1
        figure(fignum+1)
        semilogy(phoutv, Htransport, 'b')
        hold on
        Htransportsalt = res.Hin_pm;
    else
        semilogy(phoutv, Cconversion, 'r')
        hold on
        Cconversionsalt = p.alpha*res.c_cyto_uM*1e3*4*pi*p.Rb^2;
    end
    
    
    
end


plot(8, Cconversionsalt, 'or')
plot(8, Htransportsalt, 'ob')
xlabel('external pH')
ylabel('CO_2 and HCO_3^- fluxes [picomoles/s]')
hleg = legend('HCO_3^- active transport flux ','CO_2 to HCO_3^- conversion flux', 'Location', 'Best');
legend('boxoff')



