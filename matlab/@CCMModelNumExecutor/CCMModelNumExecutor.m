classdef CCMModelNumExecutor
    % class that executes the CCM model with given params
    % stores and outputs data appropriately.
    properties
        ccm_params;     % CCMParams instance for running the model.
        xnum =200; % number of spatial discritizations
        finaltime = 1000;        % total simulation time
        abstol = 1e-12;          % error tolerance
        reltol = 1e-13;
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
        function results = RunNumericalin(obj)

            %==========================================================================
            % Non-linear Diffusion Equation Solver
            % uses matlab ode functions
            %
            % concentrations returned are non-dimensional and need to be converted.
            %==========================================================================
            tic;
     
            %obj.ccm_params.pHoff = 1; % turn pH dependence of RuBisCO off to reduce computation time
            
            %==========================================================================
            % call ode solver obj contains all model parameters needed by
            % spherediffssnondim
            %==========================================================================
            
            options = odeset('RelTol',obj.reltol, 'AbsTol', obj.abstol);

            [t, u]=ode15s(@spherediffssnondim, obj.time, obj.u0, options, obj);
            [fintime, junk] = size(t);
            % check if reached steady state, if not run the same again
            change = abs(u(end,end)-u(end-1, end))./u(end, end);
            change2 = abs(u(end, end-1) - u(end-1, end-1))./u(end, end-1);
            while (change > 1e-5) && (change2 > 1e-5)
                disp('did not run long enough, running again')
                disp('could change the finaltime in CCMModelNumExecutor')
                [t, u]=ode15s(@spherediffssnondim, obj.time, u(end,:), options, obj);
                [fintime, junk] = size(t);
                change = abs(u(end,end)-u(end-1, end))./u(end, end);
                change2 = abs(u(end, end-1) - u(end-1, end-1))./u(end, end-1);
            end
                
            %==========================================================================
            % Change variables back
            %==========================================================================
            
            [r1,h1,c1] = variablechangesep(obj,u);
            
            results.c_mM = c1 * obj.ccm_params.Kca * 1e-3;
            results.h_mM = h1 * obj.ccm_params.Kba * 1e-3;
            results.c_csome_mM = results.c_mM(end,1);
            results.c_csome_uM = results.c_mM(end,1)*1e3;
            results.h_csome_mM = results.h_mM(end,1);
            results.h_csome_uM = results.h_mM(end,1)*1e3;
            results.h_nondim = h1;
            results.c_nondim = c1;
            results.ccm_params = obj.ccm_params;
            results.r = r1;
            results.t = t;
            results.fintime = fintime;
            results.x = obj.x;
            results.dx = obj.dx;
            results.xnum = obj.xnum;
            results.u0 = obj.u0;
         
            toc;
            
        end
        % =================================================================
        % function to parse output of ode solver
        % =================================================================
        
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
        % Make the discretized grid to solve equations on.
        % For solving the inner sphere only.
        % Grid spacing is for flux, concentration sits in between grid spacings.
        % This is why there are xnum+1 grid spacings and xnum concentration pts.
        %==========================================================================
        
        function result = get.x(obj)
            p = obj.ccm_params;
            r1 = 1; % we nondimensionalized by Rc (or Rb for the no Csome case)
            
            %switch to m coordinates
            m1 = (r1^3)/3;
            
            % make grid
            result = linspace(0, m1, obj.xnum);
        end
        %==========================================================================
        % calculate the grid spacing for the spherical coordinates
        %==========================================================================
        
        function result = get.dx(obj)
            p = obj.ccm_params;
            r1 = 1; % we nondimensionalized by Rc (or Rb for the no Csome case)
            
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
             result = linspace(0,obj.finaltime,10);
        end

    end
    
    
end

