
% Generate active HCO3- vs carboxysome permeability behavior space at
% internal pH = 8.

p = CCMParams_Csome;

kvec = logspace(-9,2, 200);
Hmaxvec = [5 10 15 20 25 30]*1000; % varying pool uM 
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
                
% change Km to Km from Whitehead paper for Cyanobium (Km  = 169.0 measured
% at pH 8
p.pHmeas = 8;
p.Km_meas = 169;

Ccrit = p.Km;
critjcRub2(ii)= (res.M.*Ccrit + p.Vmax*Ccrit.*res.P*p.Rc^3./(3*p.D*(Ccrit+p.Km))- ...
            p.kmC*p.Cout*(p.kmH_in*p.GH +p.alpha*p.GC +p.D/p.Rb^2))./...
           (p.Hout*((p.kmC+p.alpha)*p.GC + p.D/p.Rb^2)) - p.kmH_out;
       
% calculate jc where Hcyto vector is defined above
p.pHmeas = 7.8;
p.Km_meas = 350;
    for jj = 1:length(Hmaxvec)
    p.kcC;
    Hmax = Hmaxvec(jj);
    jc_Hmax(ii, jj) = p.CalcOptimalJc(Hmax);
    end
        
    

    
end



figure(116)
loglog(critjc, kvec, '-k')
hold on
loglog(critjcRub, kvec, '-r')
loglog(critjcRub2, kvec, '-c')
loglog(jc_Hmax, kvec, '--k')
xlabel('Active HCO_3^- transport, j_c, cm/s)')
ylabel('Carboxysome permeability')
legend('CA saturates','RuBisCO saturates, Km =276 @ pH=8', ...
        'RuBisCO saturates, Km = 169 @ pH =8', ...
        'cytosolic HCO_3^- pool = 5 mM', ...
        'cytosolic HCO_3^- pool = 10 mM', ...
        'cytosolic HCO_3^- pool = 15 mM', ...    
        'cytosolic HCO_3^- pool = 20 mM', ...  
        'cytosolic HCO_3^- pool = 25 mM', ...   
        'cytosolic HCO_3^- pool = 30 mM')

