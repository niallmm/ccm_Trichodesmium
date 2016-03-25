function jcSrec = CleakCsomezero(p,jc,Srec)
    p.jc = jc;
    exec = FullCCMModelExecutor(p);
    res = exec.RunAnalytical();
    jcSrec  = res.Ccsomeleak_pm/res.Cleak_pm-Srec;
end