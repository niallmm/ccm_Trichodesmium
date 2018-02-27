% Check Gradients across cell
    p = CCMParams_Csome;
 %  p.jc = 1e-4;
   p.jc = 0; 
  % p.alpha= p.jc;
  p.alpha= 1e-3;
    p.kcC = 1e-4;
    p.kcH = p.kcC;
    p.pH = 8;
%     p.alpha =0;

%     Hcytop = @(jc) calcHcytoDiff_Csome(jc, p, Hmax);
%     p.jc = fzero(Hcytop, 1e-2);
    
    exec = FullCCMModelExecutor(p);
    res = exec.RunAnalytical();
    
    Ccsome = res.c_csome_mM;
    Hcsome = res.h_csome_mM;
    
    rb = linspace(p.Rc, p.Rb, 1e3);
    

Ccyto = (p.kmC*p.Cout - (p.alpha+p.kmC)*res.c_csome_uM)*(p.D/(p.kcC*p.Rc^2)+1/p.Rc -1./rb)/...
    ((p.alpha+p.kmC)*p.GC + p.D/p.Rb^2) + res.c_csome_uM;
Hcyto = ((p.jc+p.kmH_out)*p.Hout + p.alpha*res.c_cyto_uM - p.kmH_in*res.h_csome_uM)*...
    (p.D/(p.kcC*p.Rc^2)+1/p.Rc -1./rb)/(p.kmH_in*p.GH + p.D/p.Rb^2) + res.h_csome_uM;
    

figure(5)
semilogy(rb, Ccyto*1e-3, 'r')
hold on
plot(rb, Hcyto*1e-3, 'b')
line([0 p.Rc],[Ccsome Ccsome], 'Color','r') 
line([0 p.Rc], [Hcsome Hcsome], 'Color', 'b')
xlabel('Cell radius (cm)')
ylabel('Concentration (mM)')
 line([p.Rc p.Rc], [1e-6 1e6], 'Color', [0.5 0.5 0.5], 'LineStyle', '-.', 'LineWidth', 3)
 xlim([0 5e-5])
 ylim([1e-5 100])
set(gca,'XTick',[0 1e-5 2e-5 3e-5 4e-5 5e-5])
set(gca,'YTick',[1e-5  1e-3  1e-1  10])
box off
