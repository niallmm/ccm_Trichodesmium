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
        
        % Document values and units here.
        VO;
        CratewO;
        OratewC;
        Hin;
        Hleak;
        Cleak;
        error;
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
            obj.CratewO = p.Vmax*C./(C+p.Km*(1+p.O/p.KO))*1e3*4*pi*p.Rc^3/3; % 1e3 converts from uM*cm^3/s to pmoles/s
            obj.OratewC = obj.VO*p.O./(p.O+p.KO*(1+C/p.Km))*1e3*4*pi*p.Rc^3/3;
            obj.Hin = p.jc*p.Hout*1e3*4*pi*p.Rb^2;  % 1e3 converts from uM*cm^3/s to pmoles/s
            obj.Hleak = p.kmH*(p.Hout - h_cyto_uM)*1e3*4*pi*p.Rb^2;
            obj.Cleak = p.kmC*(p.Cout - c_cyto_uM)*1e3*4*pi*p.Rb^2;
            obj.error = obj.OratewC/(obj.CratewO + obj.OratewC);
        end
        
    end
    
end

