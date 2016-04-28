
% Generate active HCO3- vs carboxysome permeability behavior space at
% internal pH = 8.

addpath(fileparts(pwd))
p = CCMParams_Csome;

% vector of carboxysome permeabilities we will solve for critical jc
kvec = logspace(-9,2, 200);
ratio = 1;

% vector of alpha values we will solve for critical jc
alphavec = [1e-6 1e-5 5e-5 8e-5 1e-4];

for jj= 1:length(alphavec)
    p.alpha = alphavec(jj);
for ii = 1:length(kvec)
    p.kcC = kvec(ii);
    p.kcH = ratio*p.kcC;
exec = FullCCMModelExecutor(p);
res = exec.RunAnalytical();

% Critical jc where CA 'becomes saturated'        

Ccrit = res.CCAsat0;
critjc(ii, jj)= (res.M.*Ccrit + p.Vmax*Ccrit.*res.P*p.Rc^3./(3*p.D*(Ccrit+p.Km))- ...
            p.kmC*p.Cout*(p.kmH_in*p.GH +p.alpha*p.GC +p.D/p.Rb^2))./...
           (p.Hout*((p.kmC+p.alpha)*p.GC + p.D/p.Rb^2)) - p.kmH_out;
        
% critical jc were Rubisco is saturated
Ccrit = p.Km;
critjcRub(ii, jj)= (res.M.*Ccrit + p.Vmax*Ccrit.*res.P*p.Rc^3./(3*p.D*(Ccrit+p.Km))- ...
            p.kmC*p.Cout*(p.kmH_in*p.GH +p.alpha*p.GC +p.D/p.Rb^2))./...
           (p.Hout*((p.kmC+p.alpha)*p.GC + p.D/p.Rb^2)) - p.kmH_out;
 
        
% calculate jc where Hcyto is defined in external script
    Hmax = 30000; %uM
    jc_Hmax(ii, jj) = p.CalcOptimalJc(Hmax);

    
    % calculate transition from CO2 recapture to CO2 facilitated uptake
    jc0 = jc_Hmax(ii,jj);
    critjc_uptake(ii, jj) = fzero(@(jc)Cleakzero(p,jc), jc0);
 
    
end
end

alphavec = [1e-6 1e-5 5e-5 8e-5 1e-4];

figure(116)
loglog(critjc, kvec, '-k')
hold on
loglog(critjcRub, kvec, '-r')
loglog(jc_Hmax, kvec, '--k')
loglog(critjc_uptake, kvec, '--m')
xlabel('Active HCO_3^- transport, j_c, cm/s)')
ylabel('Carboxysome permeability')
legend('\alpha = 10^{-6}, CA saturates','\alpha = 10^{-5}, CA saturates',...
    '\alpha = 5X10^{-5}, CA saturates', '\alpha = 8X10^{-5}, CA saturates', ...
    '\alpha = 10^{-4}, CA saturates', ...
    '\alpha = 10^{-6}, RuBisCO saturates','\alpha = 10^{-5}, RuBisCO saturates',...
    '\alpha = 5X10^{-5}, RuBisCO saturates', '\alpha = 8X10^{-5}, RuBisCO saturates', ...
    '\alpha = 10^{-4}, RuBisCO saturates',...
    '\alpha = 10^{-6}, cytosolic HCO_3^- = 30mM ','\alpha = 10^{-5}, cytosolic HCO_3^- = 30mM ',...
    '\alpha = 5X10^{-5}, cytosolic HCO_3^- = 30mM ', '\alpha = 8X10^{-5},  cytosolic HCO_3^- = 30mM ', ...
    '\alpha = 10^{-4},  cytosolic HCO_3^- = 30mM ',...
    '\alpha = 10^{-6}, faciliatated uptake =0 ','\alpha = 10^{-5}, faciliatated uptake =0 ',...
    '\alpha = 5X10^{-5}, faciliatated uptake =0 ', '\alpha = 8X10^{-5},  faciliatated uptake =0 ', ...
    '\alpha = 10^{-4},  faciliatated uptake =0')
