
%  calculate the lines of constant concentration for varying pH in the
%  carboxysome and transport

p = CCMParams_Csome;
% jc at optimal value for kc = 3e-5;
jc_opt =2.0992e-5;
p.jc = jc_opt;
exec = FullCCMModelExecutor(p);
res = exec.RunAnalytical();
CratewO = res.CratewO_um;
lowCrate = 0.1*CratewO; % a tenth of baseline rate
        
pHcsome = linspace(5, 9, 50);
Hmaxvec = [5, 10, 20, 30]*1000;


for ii = 1:length(pHcsome)
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
       
p.jc = critjcRub(ii);
exec = FullCCMModelExecutor(p);
res = exec.RunAnalytical();

costpertransport = 4;
% Total cost of the system in H+/fixation for each case.
total_Rub(ii) = totalProtonCost_CCM(res.Hin_pm, res.CratewO_pm, res.OratewC_pm, ...
                    res.OHrate_pm, costpertransport);
                
% change Km to Km from Whitehead paper for Cyanobium (Km  = 169.0 measured
% at pH 8
p.pHmeas = 8;
p.Km_meas = 169;

Ccrit = p.Km;
critjcRub2(ii)= (res.M.*Ccrit + p.Vmax*Ccrit.*res.P*p.Rc^3./(3*p.D*(Ccrit+p.Km))- ...
            p.kmC*p.Cout*(p.kmH_in*p.GH +p.alpha*p.GC +p.D/p.Rb^2))./...
           (p.Hout*((p.kmC+p.alpha)*p.GC + p.D/p.Rb^2)) - p.kmH_out;
       
p.jc = critjcRub2(ii);
exec = FullCCMModelExecutor(p);
res = exec.RunAnalytical();

costpertransport = 4;
% Total cost of the system in H+/fixation for each case.
total_Rub2(ii) = totalProtonCost_CCM(res.Hin_pm, res.CratewO_pm, res.OratewC_pm, ...
                    res.OHrate_pm, costpertransport);
                

p.pHmeas = 7.8;
p.Km_meas = 350;
% critical jc where CO2 fixation is 10% of baseline rate (we assume the
% RuBisCO is unsaturated.. a reasonable assumption.
Ccrit = lowCrate*p.Km*(1+p.O/p.KO)/(p.Vmax*p.Vcsome*1e-3);
critjclowC(ii)= (res.M.*Ccrit + p.Vmax*Ccrit.*res.P*p.Rc^3./(3*p.D*(Ccrit+p.Km))- ...
            p.kmC*p.Cout*(p.kmH_in*p.GH +p.alpha*p.GC +p.D/p.Rb^2))./...
           (p.Hout*((p.kmC+p.alpha)*p.GC + p.D/p.Rb^2)) - p.kmH_out;       
p.jc = critjclowC(ii);
exec = FullCCMModelExecutor(p);
res = exec.RunAnalytical();


costpertransport = 4;
% Total cost of the system in H+/fixation for each case.
total_lowC(ii) = totalProtonCost_CCM(res.Hin_pm, res.CratewO_pm, res.OratewC_pm, ...
                    res.OHrate_pm, costpertransport);

% calculate jc where Hcyto is defined in external script
%     Hmax = 30000; %uM
    for jj = 1:length(Hmaxvec)
    p.kcC;
    Hmax = Hmaxvec(jj);
    jc_Hmax(ii, jj) = p.CalcOptimalJc(Hmax);
    p.jc = jc_Hmax(ii, jj);
    exec = FullCCMModelExecutor(p);
    res = exec.RunAnalytical();
    
costpertransport = 4;
% Total cost of the system in H+/fixation for each case.
total_Hmax(ii,jj) = totalProtonCost_CCM(res.Hin_pm, res.CratewO_pm, res.OratewC_pm, ...
                    res.OHrate_pm, costpertransport);

    end
    
end


tag = {'--','-', ':', '-.'};
figure(21)
% loglog(pHcsome,critjc, '-k')
loglog(pHcsome, critjcRub, '-r')
hold on
loglog(pHcsome, critjcRub2, '-m')

loglog(pHcsome, critjclowC, '-c')
for ii = 1:length(Hmaxvec)
loglog(pHcsome, jc_Hmax(:,ii),[tag{ii} 'k'])
end
ylabel('Active HCO_3^- transport, j_c, cm/s)')
xlabel('Carboxysome pH')
legend('RuBisCO saturated',...
    'RuBisCO saturated K_m = 169 mM', ...
    'CO_2 fixation at 10%',...
    '5 mM Ci pool',...
    '10 mM Ci pool',...
    '20 mM Ci pool',...    
    '30 mM Ci pool')
legend('boxoff')

figure(121)
loglog(pHcsome, total_Rub, 'r')
hold on
loglog(pHcsome, total_lowC, 'c')
loglog(pHcsome, total_Hmax, '--k')
xlabel('Carboxysome pH')
ylabel('Cost of CCM [H^+ per CO_2 fixation]')
