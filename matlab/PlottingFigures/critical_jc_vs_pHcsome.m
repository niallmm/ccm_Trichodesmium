
%  calculate the lines of constant concentration for varying pH in the
%  carboxysome and transport

p = CCMParams_Csome;
pHcsome = linspace(5, 9, 50);
Hmaxvec = [5, 10, 15, 20, 25, 30]*1000;

  ratio = 1;


for ii = 1:length(pHcsome)
%     p.kcC = kvec(ii);
%     p.kcH = ratio*p.kcC;
p.pH_csome = pHcsome(ii);
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
%     Hmax = 30000; %uM
    for jj = 1:length(Hmaxvec)
    p.kcC;
    Hmax = Hmaxvec(jj);
    jc_Hmax(ii, jj) = p.CalcOptimalJc(Hmax);
    end
    
end



figure(21)
loglog(critjc, pHcsome, '-k')
hold on
loglog(critjcRub, pHcsome, '-r')
loglog(jc_Hmax, pHcsome, '--k')
xlabel('Active HCO_3^- transport, j_c, cm/s)')
ylabel('Carboxysome pH')
