ksweep = logspace(-9, 2, 1e3);
addpath('/Users/niallmangan/GitHub/ccm/matlab')
p = CCMParams_Csome;
p.Hout = 15;
p.pH = 8;
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

figure(5)
hold on
loglog(jc*p.Hout*4*pi*p.Rb^2*1e6, ksweep, '--k')
