
% define p = CCMParams_Csome; in external script.


% exec = FullCCMModelExecutor(p);
% res = exec.RunAnalytical();
%        

p = CCMParams_Csome;
% p.Cout = 200;
% p.NRub= 1000;
p.Km7_8 = 169;
kvec = logspace(-9,1, 200);
Hmaxvec = [5, 10, 15, 20, 25, 30]*1000;
%  p.alpha = 5e-5;
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
        
% calculate jc where Hcyto is defined in external script
%     Hmax = 30000; %uM
    for jj = 1:length(Hmaxvec)
    p.kcC;
    Hmax = Hmaxvec(jj);
    jc_Hmax(ii, jj) = p.CalcOptimalJc(Hmax);
    end
    
%     % calculate transition from CO2 recapture to CO2 facilitated uptake
%     jc0 = jc_Hmax(ii,1);
%     critjc_uptake(ii) = fzero(@(jc)Cleakzero(p,jc), jc0);
    
%     %calculate the transition where some level of CO2 is recovered by
%     %scavenging
%     Srec = 0.9;% percent of CO2 recovered by scavenging
%     critjc_Srec(ii)= fzero(@(jc)CleakCsomezero(p,jc,Srec), jc0);
    
    
%     Ccrit = p.Cout*((p.alpha+p.kmC)*p.GC*p.Rb^2/p.D + 1 -p.kmC*p.Rb^2/p.D);
%     
%     critjc_uptake(ii)= (res.M.*Ccrit + p.Vmax*Ccrit.*res.P*p.Rc^3./(3*p.D*(Ccrit+p.Km))- ...
%             p.kmC*p.Cout*(p.kmH_in*p.GH +p.alpha*p.GC +p.D/p.Rb^2))./...
%            (p.Hout*((p.kmC+p.alpha)*p.GC + p.D/p.Rb^2)) - p.kmH_out;
%     % calculate CO2 level at which there is 50% scavenging efficiency
%     Seff = 0.50;
%     Ccrit = p.Cout*(p.kmC*Seff*((p.alpha+p.kmC)*p.GC + p.D/p.Rb^2)...
%             - p.kmC*(p.kmC*Seff*p.GC + p.D/p.Rc^2))./ ...
%             (p.kmC*Seff*((p.alpha+p.kmC)*p.GC+p.D/p.Rb^2)  ...
%             -(p.alpha+p.kmC)*(p.kmC*Seff*p.GC + p.D/p.Rc^2));
%     critjc_50S(ii) = (res.M.*Ccrit + p.Vmax*Ccrit.*res.P*p.Rc^3./(3*p.D*(Ccrit+p.Km))- ...
%             p.kmC*p.Cout*(p.kmH_in*p.GH +p.alpha*p.GC +p.D/p.Rb^2))./...
%            (p.Hout*((p.kmC+p.alpha)*p.GC + p.D/p.Rb^2)) - p.kmH_out;
       
            

        
    

    
end


% figure(6)        
% loglog(critjc*p.Hout*4*pi*p.Rb^2*1e3, kvec, '-k')
% hold on
% loglog(critjcRub*p.Hout*4*pi*p.Rb^2*1e3, kvec, '-r')
% %plot(critjcRub99*p.Hout*4*pi*p.Rb^2*1e6, kvec, '-c')
% loglog(jc_Hmax*p.Hout*4*pi*p.Rb^2*1e3, kvec, '--k')
% xlabel('Active HCO_3^- transport, j_c*H_{out}, (picomole/s)')
% ylabel('Carboxysome permeability')
% 
% figure(116)
% loglog(critjc, kvec, '-k')
% hold on
% loglog(critjcRub, kvec, '-r')
% loglog(jc_Hmax, kvec, '--k')
% %plot(critjcRub99, kvec, '-c')
% 
% xlabel('Active HCO_3^- transport, j_c, cm/s)')
% ylabel('Carboxysome permeability')
