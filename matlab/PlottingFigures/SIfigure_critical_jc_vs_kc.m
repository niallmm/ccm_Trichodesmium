
% Generate active HCO3- vs carboxysome permeability behavior space at
% internal pH = 8.

p = CCMParams_Csome;

kvec = logspace(-9,2, 200);
ratio = 1;


for ii = 1:length(kvec)
    p.kcC = kvec(ii);
    p.kcH = ratio*p.kcC;
exec = FullCCMModelExecutor(p);
res = exec.RunAnalytical();

% Critical jc where CA 'becomes saturated'        

Ccrit = res.CCAsat0;
critjc(ii)= (res.M.*Ccrit + p.Vmax*Ccrit.*res.P*p.Rc^3./(3*p.D*(Ccrit+p.Km))- ...
            p.kmC*p.Cout*(p.kmH_in*p.GH +p.alpha*p.GC +p.D/p.Rb^2))./...
           (p.Hout*((p.kmC+p.alpha)*p.GC + p.D/p.Rb^2)) - p.kmH_out;
        
% critical jc were Rubisco is saturated
Ccrit = p.Km;
critjcRub(ii)= (res.M.*Ccrit + p.Vmax*Ccrit.*res.P*p.Rc^3./(3*p.D*(Ccrit+p.Km))- ...
            p.kmC*p.Cout*(p.kmH_in*p.GH +p.alpha*p.GC +p.D/p.Rb^2))./...
           (p.Hout*((p.kmC+p.alpha)*p.GC + p.D/p.Rb^2)) - p.kmH_out;
 
        
% calculate jc where Hcyto is defined in external script
     Hmax = 30000; %uM
    Hmax = Hmaxvec(jj);
    jc_Hmax(ii, jj) = p.CalcOptimalJc(Hmax);

       
    

    
end

grey = [0.4, 0.4, 0.4];

% 
figure(116)
loglog(critjc, kvec, '-k')
hold on
loglog(critjcRub, kvec, '-r')
loglog(jc_Hmax, kvec, 'Color', grey, 'LineStyle', '--')
line([1e-6 1e3],[ 3e-5 3e-5], 'Color', 'k', 'LineStyle', '--')
xlabel('Active HCO_3^- transport, j_c, cm/s)')
ylabel('Carboxysome permeability')
legend('CA saturated','RuBisCO saturated',... 
    '30 mM Ci pool')
legend('boxoff')