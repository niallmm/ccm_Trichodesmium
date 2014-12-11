classdef FullCCMAnalyticalSolution
    % Concentrations are in (TxN) matrices where T is the number of
    % timepoints and N is the number of points along the radius of the cell
    % considered in the discretization of the cell.
    properties
        ccm_params;     % params used to solve the model
        h_cyto_uM;      % uM concentration of total bicarbonate in cytoplasm.
        c_cyto_uM;      % uM concentration of CO2 in cytoplasm.
        h_cyto_mM;      % mM concentration of total bicarbonate in cytoplasm.
        c_cyto_mM;      % mM concentration of CO2 in cytoplasm.
        h_csome_uM;     % uM concentration of total bicarbonate in carboxysome.
        c_csome_uM;     % uM concentration of CO2 in carboxysome.
        h_csome_mM;     % mM concentration of total bicarbonate in carboxysome.
        c_csome_mM;     % mM concentration of CO2 in carboxysome.
        
        % =================================================================
        % Calculate CO2 and O2 fixation rates for whole cell
        VO;             % [uM/s] maximum rate of oxygen fixation calculated from the specificity of RuBisCO 
        % intgrated over carboxysome volume 
        CratewO_um;        % [um/s] rate of CO2 fixation with oxygen accounted for
        OratewC_um;        % [um/s] rate of O2 fication with CO2 accounted for
        CratewO_pm;        % [pmole/s] rate of CO2 fixation with oxygen accounted for
        OratewC_pm;        % [pmole/s] rate of O2 fication with CO2 accounted for
        % =================================================================
        % Calculate CO2 and HCO3- flux rates at cell membrane
        % integrated over surface area of cell 
        Hin_pm;            % [pmole/s] rate of active uptake of HCO3- jc*Hout
        Hleak_pm;          % [pmole/s] rate of HCO3- leakage out of cell kmH*(Hout-Hcyto)
        Cleak_pm;          % [pmole/s] rate of CO2 leakage out of cell kmC*(Cout-Ccyto)
        Hin_um;            % [um/s] 
        Hleak_um;          % [um/s]
        Cleak_um;          % [um/s]
        
        error;             % the proportion of oxygen fixations to total fixation events
    end
    
    methods
        function obj = FullCCMAnalyticalSolution(ccm_params, h_cyto_uM, c_cyto_uM, h_csome_uM, c_csome_uM)
            obj.ccm_params = ccm_params;
            obj.h_cyto_uM = h_cyto_uM;
            obj.c_cyto_uM = c_cyto_uM;
            obj.h_cyto_mM = h_cyto_uM * 1e-3;
            obj.c_cyto_mM = c_cyto_uM * 1e-3;
            obj.h_csome_uM = h_csome_uM;
            obj.c_csome_uM = c_csome_uM;
            obj.h_csome_mM = h_csome_uM * 1e-3;
            obj.c_csome_mM = c_csome_uM * 1e-3;
            
            p = ccm_params;
            C = obj.c_csome_uM;
            obj.VO = p.Vmax*p.KO/(p.Km*p.S_sat);
            obj.CratewO_pm = p.Vmax*C./(C+p.Km*(1+p.O/p.KO))*p.Vcsome*1e3; 
            obj.CratewO_um = p.Vmax*C./(C+p.Km*(1+p.O/p.KO))*p.Vcsome*1e-3; % convert from uM*cm^3 to umoles
            obj.OratewC_pm = obj.VO*p.O./(p.O+p.KO*(1+C/p.Km))*p.Vcsome*1e3;
            obj.OratewC_um = obj.VO*p.O./(p.O+p.KO*(1+C/p.Km))*p.Vcsome*1e-3;
            
            obj.Hin_pm = p.jc*p.Hout*p.SAcell*1e3; 
            obj.Hleak_pm = p.kmH*(p.Hout - h_cyto_uM)*p.SAcell*1e3;
            obj.Cleak_pm = p.kmC*(p.Cout - c_cyto_uM)*p.SAcell*1e3;
            obj.Hin_um = p.jc*p.Hout*p.SAcell*1e-3; 
            obj.Hleak_um = p.kmH*(p.Hout - h_cyto_uM)*p.SAcell*1e-3;
            obj.Cleak_um = p.kmC*(p.Cout - c_cyto_uM)*p.SAcell*1e-3;
            
            obj.error = obj.OratewC_pm/(obj.CratewO_pm + obj.OratewC_pm);
        end
        
    end
    
end

