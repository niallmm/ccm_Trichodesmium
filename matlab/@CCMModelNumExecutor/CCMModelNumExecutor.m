classdef CCMModelNumExecutor
    % class that executes the CCM model with given params
    % stores and outputs data appropriately.
    properties
        ccm_params;     % CCMParams instance for running the model.
        xnum = 50; % number of spatial discritizations
        finaltime = 1000;        % total simulation time
        abstol = 1e-5;          % error tolerance
        reltol = 1e-8;
    end
    
    properties (Dependent)
        u0 % intial conditions
        x   % initialize the discritized x vector (set in CCMModelNumExecutor)
        dx  % initialize the discritized spacing. (set in CCMModelNumExecutor)
        time
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
            %
            % concentrations returned are non-dimensional and need to be converted.
            %==========================================================================
            tic;
            %==========================================================================
            % set parameters (contains model parameters, location of enzymes
            % (in compartment or throughout cell)
            %==========================================================================
                  
            obj.ccm_params.pHoff = 1; % turn pH dependence of RuBisCO off to reduce computation time
            
            %==========================================================================
            % call ode solver
            %==========================================================================
            
            options = odeset('RelTol',obj.reltol, 'AbsTol', obj.abstol);

            [t, u]=ode15s(@spherediffssnondim, obj.time, obj.u0, options, obj);
            [fintime, junk] = size(t);
            %==========================================================================
            % Change variables back
            %==========================================================================
            
            [r1,h1,c1] = variablechangesep(obj,u);
            
            results.c_mM = c1 * obj.ccm_params.Kca * 1e-3;
            results.h_mM = h1 * obj.ccm_params.Kba * 1e-3;
            results.h_nondim = h1;
            results.c_nondim = c1;
            results.ccm_params = obj.ccm_params;
            results.r = r1;
            results.t = t;
            results.fintime = fintime;
            toc;
            
        end
        
        [r1,h1,c1] = variablechangesep(obj,u);

        %==========================================================================
        % call to set an initial condition only in the carboxysome
        %==========================================================================
        
        function result = get.u0(obj)
            initv = zeros(2, obj.xnum); % right now the initial condition is just set ot zero, but we could change that.
            h0 = initv(1,:);
            c0 = initv(2,:);
            
            result(1:2:obj.xnum*2)= h0;
            
            result(2:2:obj.xnum*2)= c0;
        end
        %==========================================================================
        % create grid in spherical radial coordinates
        %==========================================================================
        
        function result = get.x(obj)
            p = obj.ccm_params;
            r1 = p.Rc;
            
            %switch to m coordinates
            m1 = (r1^3)/3;
            
            % grid spacing
            dx1 = m1/(obj.xnum);
            
            % make grid
            result = linspace(0, m1, obj.xnum);
        end
        %==========================================================================
        % calculate the grid spacing for the spherical coordinates
        %==========================================================================
        
        function result = get.dx(obj)
            p = obj.ccm_params;
            r1 = p.Rc;
            
            %switch to m coordinates
            m1 = (r1^3)/3;
            
            % grid spacing
            result = m1/(obj.xnum);
        end
        %==========================================================================
        % create a vector of time points at which to request the solution
        % from ode. Generally we are just using ode to make sure we are
        % at steady state, not to find time dependent solutions, but it is
        % possible. ode will calculate at whatever time points needed to
        % meet the error tolerances, but will only save them for the values
        % in this vector.
        %==========================================================================
        
        function result = get.time(obj)
             result = linspace(0,obj.finaltime,100);
        end

    end
    
    
end

