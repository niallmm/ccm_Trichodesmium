% Soluions w/ reactions in whole cell
% Concentrations returned are in uM units.
function [H, C] = CHconc_NoCsome(jc_in, ccm_params)

p = ccm_params;


% jc = jcsweep;
% CA equilibriates + Rub
N2 = p.kmC*p.Cout+(p.jc+p.kmH_out)*p.Hout;
M2 = p.kmC*(1+p.kmH_in/p.Keq);

CcytoRub = 0.5*(N2/M2 - p.Rb*p.Vmax/(3*M2) - p.Km)...
          +0.5*sqrt((p.Km-N2/M2 + p.Rb*p.Vmax/(3*M2)).^2 + 4*p.Km*N2/M2);

HcytoRub= CcytoRub/p.Keq;
      
CcytoRubsat = N2/M2 - p.Rb*p.Vmax/3;

CcytoRubunsat = N2/(M2+p.Rb*p.Vmax/(3*p.Km));

% CA saturated

CcytoCAsat0 = p.kmC*p.Cout/(p.alpha+p.kmC) + p.Vba*(p.Rb/(3*(p.alpha+p.kmC))+ p.Rb^2/(6*p.D));

CcytoCAsatRb = p.kmC*p.Cout/(p.alpha+p.kmC) + p.Vba*(p.Rb/(3*(p.alpha+p.kmC)));

HcytoCAsat0 = -p.Vba*(p.Rb^2/(6*p.D) + p.Rb/(3*p.kmH_in)) + (p.jc+p.kmH_out)*p.Hout/p.kmH_in ...
                + p.alpha*CcytoCAsatRb/p.kmH_in;

            
if  CcytoRub>=CcytoCAsat0
    C = CcytoCAsat0;
    H = HcytoCAsat0;
    %csat = 1
elseif CcytoRub<CcytoCAsat0
    C = CcytoRub;
    H = HcytoRub;
    %CAunsat =1
end

