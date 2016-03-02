classdef NoCsomeModelExecutor 
        properties
        ccm_params;     % CCMParams instance for running the model.
    end
    
    methods
        function obj = NoCsomeModelExecutor(ccm_params)
              obj.ccm_params = ccm_params;
        end
        
        function result = RunAnalytical(obj)
            p = obj.ccm_params;
            result = NoCsomeAnalyticalSolution(obj.ccm_params);
        end
        function  resultN = RunNumerical(obj)
             exec= CCMModelNumExecutor(obj.ccm_params);
            resultN = exec.RunNumericalin();
        end
    end
    
end

