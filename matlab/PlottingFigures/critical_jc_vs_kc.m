
addpath('/Users/niallmangan/GitHub/ccm/matlab')
p = CCMParams_Csome;
p.pH = 8;


exec = FullCCMModelExecutor(p);
res = exec.RunAnalytical();
       
        
kvec = logspace(-9,2, 1e3);

for i = 1:length(kvec)
    p.kcC = kvec(i);
    p.kcH = p.kcC;
exec = FullCCMModelExecutor(p);
res = exec.RunAnalytical();



% Critical jc where CA 'becomes saturated'        

Ccrit = res.CCAsat0;
critjc(i)= (res.M.*Ccrit + p.Vmax*Ccrit.*res.P*p.Rc^3./(3*p.D*(Ccrit+p.Km))- ...
            p.kmC*p.Cout*((p.kmH_in+p.alpha)*p.GH +p.D/p.Rb^2))./...
           (p.Hout*((p.kmC+p.alpha)*p.GC + p.D/p.Rb^2)) - p.kmH_out;
        
% critical jc were Rubisco is saturated
Ccrit = p.Km;
critjcRub(i)= (res.M.*Ccrit + p.Vmax*Ccrit.*res.P*p.Rc^3./(3*p.D*(Ccrit+p.Km))- ...
            p.kmC*p.Cout*((p.kmH_in+p.alpha)*p.GH +p.D/p.Rb^2))./...
           (p.Hout*((p.kmC+p.alpha)*p.GC + p.D/p.Rb^2)) - p.kmH_out;
        
% % calculate jc where carboxylation is 99% 
% % when RuBisCO is saturated: S = vc/vo ~ 13;
% % not saturated: S = (vc/Kc)/(vo/Ko)~ 50 %SavirMilo paper
% % carboxylation/oxyp.Genation = S*[CO_2]/[O_2] -> 
% % CO_2 needed to make carboxylation 99% -> 99*[O_2]/S
% O2 = p.O; % uM (Savir Milo paper)
% S = p.S_sat; % when RuBisCO is saturated
% 
% Ccrit = 99*O2/S;
% critjcRub99(i)= (res.M.*Ccrit + p.Vmax*Ccrit.*res.P*p.Rc^3./(3*p.D*(Ccrit+p.Km))- ...
%             p.kmC*p.Cout*((p.kmH_in+p.alpha)*p.GH +p.D/p.Rb^2))./...
%            (p.Hout*((p.kmC+p.alpha)*p.GC + p.D/p.Rb^2)) - p.kmH_out;
        
% calculate jc where Hcyto = 30mM
    Hmax = 30000; %uM
    p.kcC
    jc_Hmax(i) = p.CalcOptimalJc(Hmax);
end
figure(6)        
loglog(critjc*p.Hout*4*pi*p.Rb^2*1e3, kvec, '-k')
hold on
loglog(critjcRub*p.Hout*4*pi*p.Rb^2*1e3, kvec, '-r')
%plot(critjcRub99*p.Hout*4*pi*p.Rb^2*1e6, kvec, '-c')
loglog(jc_Hmax*p.Hout*4*pi*p.Rb^2*1e3, kvec, '--k')
xlabel('Active HCO_3^- transport, j_c*H_{out}, (picomole/s)')
ylabel('Carboxysome permeability')

figure(116)
loglog(critjc, kvec, '-k')
hold on
loglog(critjcRub, kvec, '-r')
loglog(jc_Hmax, kvec, '--k')
%plot(critjcRub99, kvec, '-c')

xlabel('Active HCO_3^- transport, j_c, cm/s)')
ylabel('Carboxysome permeability')
