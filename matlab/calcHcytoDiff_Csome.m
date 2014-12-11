% Calculate the difference between the cytoplasmic bicarbonate
% concentration and the maximum supplied. Used for finding optimal jc.
function Hdiff = calcHcytoDiff_Csome(jc, ccm_params, Hmax)
    [HcytoRb, CcytoRb, H, C] = CHconc_Csome(jc, ccm_params);
    Hdiff = HcytoRb-Hmax;