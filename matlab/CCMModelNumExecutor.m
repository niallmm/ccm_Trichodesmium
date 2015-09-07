classdef CCMModelNumExecutor
    % class that executes the CCM model with given params
    % stores and outputs data appropriately.
    properties
        ccm_params;     % CCMParams instance for running the model.
    end
    
    methods
        function obj = CCMModelNumExecutor(ccm_params)
            obj.ccm_params = ccm_params;
        end
        
        % Runs a numerical simulation of the CCM system to find the
        % steady state behavior. Stores results in this class.
        % Note: numerical code is generic to the various cases we consider
        % (e.g. no carboxysome, no ccm, etc) so this method can be
        % implemented generically.
        function results = RunNumerical(obj)
            %==========================================================================
            % Non-linear Diffusion Equation Solver
            % uses matlab ode functions
            % calls setgrid.m, initold.m, variablechangesep.m
            % cytosol.m calculates solutions outside csome
            % checks error
            % plots solution
            % for two concentrations
            %
            % concentrations returned are non-dimensional and need to be converted.
            %==========================================================================
            tic;
            %==========================================================================
            % set parameter (contains model parameters, location of enzymes
            % (in compartment or throughout cell, and numerical parameters)
            %==========================================================================
            
            p = obj.ccm_params;  % shorthand

            %==========================================================================
            % call ode solver
            %==========================================================================
            
            options = odeset('RelTol',p.reltol, 'AbsTol', p.abstol);
=
            [t, u]=ode15s(@spherediffssnondim, p.time, p.u0, options, p);
            [fintime, junk] = size(t);
            %==========================================================================
            % Change variables back
            %==========================================================================
            
            [r1,h1,c1] = variablechangesep(x,u);
            
            toc;
            results = NumericalCCMModelSolution(p, r, h, c, fintime, t);
        end
    end
    
    methods (Abstract)
        % Run the analytical model
        RunAnalytical(obj)
    end
    
end

