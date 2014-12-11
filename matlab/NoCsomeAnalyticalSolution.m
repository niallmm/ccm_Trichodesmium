classdef NoCsomeAnalyticalSolution
    % Concentrations are in (TxN) matrices where T is the number of
    % timepoints and N is the number of points along the radius of the cell
    % considered in the discretization of the cell.
    properties
        ccm_params;     % params used to solve the model
        h_cyto_uM;      % uM concentration of total bicarbonate.
        c_cyto_uM;      % uM concentration of total bicarbonate.
        h_cyto_mM;      % mM concentration of total bicarbonate.
        c_cyto_mM;      % mM concentration of co2.
        
        % Document values and units here.
        VO;
        CratewO;
        OratewC;
        Hin;
        Hleak;
        Cleak;
        Ccyto;
        Hcyto;
        error;
    end
    
    methods
        function obj = NoCsomeAnalyticalSolution(ccm_params, h_cyto_uM, c_cyto_uM)
            obj.ccm_params = ccm_params;
            obj.h_cyto_uM = h_cyto_uM;
            obj.c_cyto_uM = c_cyto_uM;
            obj.h_cyto_mM = h_cyto_uM * 1e-3;
            obj.c_cyto_mM = c_cyto_uM * 1e-3;
            
            % Uses analytic solution in cytosol to break down the fate of carbon in the
            % case that all the CO2 fixation maxhinery is in the cytosol.
            p = ccm_params;

            % CA saturated case
            CcytoCAsat0 = p.kmC*p.Cout/(p.alpha+p.kmC) + p.VbaCell*(p.Rb/(3*(p.alpha+p.kmC))+ p.Rb^2/(6*p.D));
            CcytoCAsatRb = p.kmC*p.Cout/(p.alpha+p.kmC) + p.VbaCell*(p.Rb/(3*(p.alpha+p.kmC)));
            HcytoCAsat0 = -p.VbaCell*(p.Rb^2/(6*p.D) + p.Rb/(3*p.kmH)) + (p.jc+p.kmH)*p.Hout/p.kmH ...
                            + p.alpha*CcytoCAsatRb/p.kmH;
            
            if  obj.c_cyto_uM >=CcytoCAsat0
                C = CcytoCAsat0;
                H = HcytoCAsat0;
                %csat = 1
            elseif obj.c_cyto_uM <CcytoCAsat0
                C = obj.c_cyto_uM;
                H = obj.h_cyto_uM;
                %CAunsat = 1
            end

            obj.VO = p.Vmax*p.KO/(p.Km*p.S_sat);
            obj.CratewO = p.Vmax*C./(C+p.Km*(1+p.O/p.KO)); % 1e3 converts from uM*cm^3/s to pmoles/s
            obj.OratewC = obj.VO*p.O./(p.O+p.KO*(1+C/p.Km));
            obj.Hin = p.jc*p.Hout*1e3*4*pi*p.Rb^2;  % 1e3 converts from uM*cm^3/s to pmoles/s
            obj.Hleak = p.kmH*(p.Hout - obj.h_cyto_uM)*1e3*4*pi*p.Rb^2;
            obj.Cleak = p.kmC*(p.Cout - obj.c_cyto_uM)*1e3*4*pi*p.Rb^2;
            obj.error = obj.OratewC/(obj.CratewO+obj.OratewC);
        end
    end
    
end

