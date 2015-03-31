% Check Gradients across cell
    p = CCMParams_Csome;
    p.jc = 3.5e-4;
    p.k = 1e-4;
    p.pH = 8;
    exec = FullCCMModelExecutor(p);
    res = exec.RunAnalytical();

    Ccsome = res.c_csome_mM;
    Hcsome = res.h_csome_mM;
    
    rb = linspace(p.Rc, p.Rb, 1e3);
    

Ccyto = (p.kmC*p.Cout - (p.alpha+p.kmC)*res.c_csome_uM)*(p.D/(p.k*p.Rc^2)+1/p.Rc -1./rb)/...
    ((p.alpha+p.kmC)*p.G + p.D/p.Rb^2) + res.c_csome_uM;
Hcyto = ((p.jc+p.kmH)*p.Hout + p.alpha*res.c_cyto_uM - p.kmH*res.h_csome_uM)*...
    (p.D/(p.k*p.Rc^2)+1/p.Rc -1./rb)/(p.kmH*p.G + p.D/p.Rb^2) + res.h_csome_uM;
    

figure
semilogy(rb, Ccyto*1e-3, 'r')
hold on
plot(rb, Hcyto*1e-3, 'b')
line([0 p.Rc],[Ccsome Ccsome], 'Color','r') 
line([0 p.Rc], [Hcsome Hcsome], 'Color', 'b')
xlabel('Cell radius (cm)')
ylabel('Concentration (mM)')
 line([Rc Rc], [1e-6 1e6], 'Color', [0.5 0.5 0.5], 'LineStyle', '-.', 'LineWidth', 3)
xlim([0 5e-5])
ylim([1e-5 100])
set(gca,'XTick',[0 1e-5 2e-5 3e-5 4e-5 5e-5])
set(gca,'YTick',[1e-5  1e-3  1e-1  10])
box off
