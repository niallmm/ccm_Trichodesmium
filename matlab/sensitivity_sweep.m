% calculate sweeps through sensitivity analysis

resultsfolder = 'savedoutput/03102016Code/';
                    mkdir(resultsfolder)
addpath('./PlottingFigures')


Rcvec = 0.5*[2, 10, 20, 40]*1e-6; % 20 nm, 100nm, 200 nm, to 400 nm diameter carboxysomes
ratiovec= [1 10 100 1000]; % ratio between kcC and kcH
kmHvec = [3e-4 3e-3 3e-2]; % membrane permeability for H2CO3

p = CCMParams_Csome;


%     for jj = 1:length(Rcvec)
%         p.Rc = Rcvec(jj);
%         for kk = 1:length(ratiovec)
%             ratio = ratiovec(kk);
               for nn = 1:length(kmHvec)
                   p.kmH_base = kmHvec(nn);


%                     
%                     critical_jc_vs_kc
                    
                    foldername = ['Rc', num2str(p.Rc), ...
                        'ratio' num2str(ratio), ...
                        'kmHbase' num2str(p.kmH_base), '/']
                    savefolder = [resultsfolder, foldername];
%                     mkdir(savefolder)
%                     save([savefolder, '/jc_vs_kc.mat'])
                    
                    load([savefolder, 'jc_vs_kc.mat'])
                    figure(116)
%                     loglog(critjc, kvec, '-k')
%                     hold on
                    loglog(critjcRub, kvec, '-r')
                    hold on
%                     loglog(jc_Hmax, kvec, '--k')
                    xlabel('Active HCO_3^- transport, j_c, cm/s)')
                    ylabel('Carboxysome permeability')
                    drawnow
                    clear exec res Ccrit critjc critjcRub jcHmax kvec

             %   end
%             end
        end


                    
                    
                    