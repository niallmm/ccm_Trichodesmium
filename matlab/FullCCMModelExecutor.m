classdef FullCCMModelExecutor %  < CCMModelExector    
    properties
        ccm_params;     % CCMParams instance for running the model.
    end
    
    methods
         function obj = FullCCMModelExecutor(ccm_params)
              obj.ccm_params = ccm_params;
         end
        
        function result = RunAnalytical(obj)
            p = obj.ccm_params;
            result = FullCCMAnalyticalSolution(p);
        end
    end
    
end

