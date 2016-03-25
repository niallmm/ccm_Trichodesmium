function jcCleak0 = Cleakzero(p,jc)
    p.jc = jc;
    exec = FullCCMModelExecutor(p);
    res = exec.RunAnalytical();
    jcCleak0  = res.Cleak_um;
end