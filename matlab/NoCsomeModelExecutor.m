classdef NoCsomeModelExecutor < CCMModelExector    
    methods
        function obj = NoCsomeModelExecutor(ccm_params)
            obj@CCMModelExector(ccm_params); 
        end
        
        function result = RunAnalytical(obj)
            p = obj.ccm_params;
            [Hcyto, Ccyto] = CHconc_NoCsome(p.jc, p);
            result = NoCsomeAnalyticalSolution(obj.ccm_params, Hcyto, Ccyto);
        end
    end
    
end

