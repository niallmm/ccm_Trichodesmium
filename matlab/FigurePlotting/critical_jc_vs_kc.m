
addpath('/Users/niallmangan/GitHub/ccm/matlab')
p = CCMParams_Csome;
p.pH = 8;
exec = FullCCMModelExecutor(p);
res = exec.RunAnalytical();
       
        
k = logspace(-9,2, 1e4);


G = D./(k*Rc^2) + 1/Rc - 1/Rb;


M = (kmC + alpha)*(1+Vca*Kba/(Vba*Kca))*kmH*G + ...
    kmC*(1+kmH*Vca*Kba/(kmC*Vba*Kca))*D/Rb^2;
P = ((alpha + kmC)*G + D/Rb^2).*(kmH*G + D/Rb^2);


CCAsat0 = Vba*(Rc^3)*(G+D/((alpha+kmC)*Rb^2))/(3*D) + ...
            Vba*(Rc^2)/(6*D) + kmC*Cout/(alpha+kmC);

% Critical jc where CA 'becomes saturated'        

Ccrit = CCAsat0;
critjc= (M.*Ccrit + Vmax*Ccrit.*P*Rc^3./(3*D*(Ccrit+Km))- ...
            kmC*Cout*((kmH+alpha)*G +D/Rb^2))./...
           (Hout*((kmC+alpha)*G + D/Rb^2)) - kmH;
        
% critical jc were Rubisco is saturated
Ccrit = Km;
critjcRub= (M*Ccrit + Vmax*Ccrit*P*Rc^3/(3*D*(Ccrit+Km))- ...
            kmC*Cout*((kmH+alpha)*G +D/Rb^2))./...
           (Hout*((kmC+alpha)*G + D/Rb^2)) - kmH;
        
% calculate jc where carboxylation is 99% 
% when RuBisCO is saturated: S = vc/vo ~ 13;
% not saturated: S = (vc/Kc)/(vo/Ko)~ 50 %SavirMilo paper
% carboxylation/oxygenation = S*[CO_2]/[O_2] -> 
% CO_2 needed to make carboxylation 99% -> 99*[O_2]/S
O2 = 260; % uM (Savir Milo paper)
S = 13; % when RuBisCO is saturated

Ccrit = 99*O2/S;
critjcRub99= (M*Ccrit + Vmax*Ccrit*P*Rc^3/(3*D*(Ccrit+Km))- ...
            kmC*Cout*((kmH+alpha)*G +D/Rb^2))./...
           (Hout*((kmC+alpha)*G + D/Rb^2)) - kmH;
        
Ccrit = 999*O2/13;
critjcRub999= (M*Ccrit + Vmax*Ccrit*P*Rc^3/(3*D*(Ccrit+Km))- ...
            kmC*Cout*((kmH+alpha)*G +D/Rb^2))./...
           (Hout*((kmC+alpha)*G + D/Rb^2)) - kmH;
figure(3)        
loglog(critjc, k, '-k')
hold on
plot(critjcRub, k, '-r')
plot(critjcRub99, k, '-c')
% plot(critjcRub999, k, '-g')
xlabel('Active HCO_3^- transport, j_c, (cm/s)')
ylabel('Carboxysome permeability')
