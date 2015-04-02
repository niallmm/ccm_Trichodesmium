% vary kcC and kcH separately and see what happens

p = CCMParams_Csome;

kcHvec = [1e-3 1e-4 1e-5];
kcCvec = logspace(-3, -5, 10);
p.pH = 8;
p.jc = 0.01;
for j = 1:length(kcHvec)
    p.kcH = kcHvec(j);
for i = 1:length(kcCvec)
    p.kcC = kcCvec(i);
    exec = FullCCMModelExecutor(p);
    res = exec.RunAnalytical();
    Keq_vec(i,j) = p.Keq;
    Ccsome(i,j) = res.c_csome_mM;
    Hcsome(i,j) = res.h_csome_mM;
end
end
% 
% figure(11)
% semilogy(pHvec, Keq_vec, 'ok')
% xlabel('pH')
% ylabel('Keq')

figure(12)
loglog(kcCvec, Ccsome, 'r')
hold on
loglog(kcCvec, Hcsome, 'b')
xlabel('kcC')
ylabel('Ccsome')