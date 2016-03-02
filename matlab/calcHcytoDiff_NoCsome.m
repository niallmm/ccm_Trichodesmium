% find critical jc of 
function Hdiff = calcHcytoDiff_NoCsome(jc_in, ccm_params, Hmax)
    ccm_params.h_cyto_exp = Hmax;
    ccm_params.jc = jc_in;
    exec = NoCsomeModelExecutor(ccm_params);
    res = exec.RunAnalytical();
    Hdiff = res.hdiff;