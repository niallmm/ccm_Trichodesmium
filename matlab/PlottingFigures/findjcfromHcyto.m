ksweep = logspace(-9, 2, 1e3);
addpath('/Users/niallmangan/GitHub/ccm/matlab')
p = CCMParams_Csome;
p.Hout = 15;
p.pH = 8;
p.kRub = 11.6; % rxns/s maximum reaction rate at single active site
p.Km_8 = 340;    % half max reaction rate of RuBisCO, uM
p.S_sat = 43;  % specificity ratio
p.KO = 972;    % uM
exec = FullCCMModelExecutor(p);
res = exec.RunAnalytical();
for i= 1:length(ksweep)
    p.k = ksweep(i);
exec = FullCCMModelExecutor(p);
res = exec.RunAnalytical();
     p.alpha =0;
    Hmax = 30000; %uM
%     Hcytop = @(jc) calcHcyto(p.jc, p.k, alpha, Hmax, kmH);(jc, ccm_params, Hmax)
    Hcytop = @(jc) calcHcytoDiff_Csome(jc, p, Hmax);
    jc(i) = fzero(Hcytop, 1e-2);
%     Hcytop(jc(i))
end

figure(6)
hold on
loglog(jc*p.Hout*4*pi*p.Rb^2*1e6, ksweep, '--k')

figure(116) 
hold on
loglog(jc, ksweep, '--k')
