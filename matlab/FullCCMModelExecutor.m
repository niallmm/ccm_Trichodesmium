classdef FullCCMModelExecutor < CCMModelExector    
    methods
        function obj = FullCCMModelExecutor(ccm_params)
            obj@CCMModelExector(ccm_params); 
        end
        
        function result = RunAnalytical(obj)
            p = obj.ccm_params;
            result = FullCCMAnalyticalSolution(p);
        end
    end
    
end

