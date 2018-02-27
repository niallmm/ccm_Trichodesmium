classdef FullCCMModelExecutor    
    properties
        ccm_params;     % CCMParams instance for running the model.
    end
    
    methods
         function obj = FullCCMModelExecutor(ccm_params)
              obj.ccm_params = ccm_params;
         end
        
        function result = RunAnalytical(obj)
            p = obj.ccm_params;
            result = FullCCMAnalyticalSolution(p);
        end
        
        function resultN = RunNumerical(obj)
            exec= CCMModelNumExecutor(obj.ccm_params);
            resultN = exec.RunNumericalin();
            p = obj.ccm_params;
 
            
            % concentration in the cytosol at r = Rb
            resultN.c_cyto_uM = (p.kmC*p.Cout - (p.alpha+p.kmC)*resultN.c_csome_uM)*p.GC/...
                ((p.alpha+p.kmC)*p.GC + p.D/p.Rb^2) +resultN.c_csome_uM;
            
            resultN.h_cyto_uM = ((p.jc+p.kmH_out)*p.Hout + p.alpha*resultN.c_cyto_uM - ...
               p.kmH_in*resultN.h_csome_uM)*p.GH/(p.kmH_in*p.GH + p.D/p.Rb^2)+resultN.h_csome_uM;
           % concentration across the cell
            resultN.r_cell = linspace(p.Rc, p.Rb, 50);
            resultN.c_cyto_rad_uM = (p.kmC*p.Cout - (p.alpha+p.kmC)*resultN.c_csome_uM)...
               *(p.D/(p.kcC*p.Rc^2) + 1/p.Rc - 1./resultN.r_cell)/...
                ((p.alpha+p.kmC)*p.GC + p.D/p.Rb^2) +resultN.c_csome_uM;
            resultN.c_cyto_rad_mM = resultN.c_cyto_rad_uM*1e-3;
            resultN.h_cyto_rad_uM = ((p.jc+p.kmH_out)*p.Hout + p.alpha*resultN.c_cyto_uM - ...
               p.kmH_in*resultN.h_csome_uM)*(p.D/(p.kcH*p.Rc^2) + 1/p.Rc - 1./resultN.r_cell)...
               /(p.kmH_in*p.GH + p.D/p.Rb^2)+resultN.h_csome_uM;
            resultN.h_cyto_rad_mM = resultN.h_cyto_rad_uM*1e-3;
        end
    end
end

