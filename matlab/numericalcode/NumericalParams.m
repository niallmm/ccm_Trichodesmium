classdef NumericalParams
    
    properties
        xnum = 100; % number of spatial discritizations
        finaltime = 100000000000;        % total simulation time
        time = linspace(0,finaltime,100);
        %time = [0 0.0001];
        abstol = 1e-12;          % error tolerance
        reltol = 1e-13;
        
    end
    
    properties (Dependent)
        u0 % intial conditions
        x   % initialize the discritized x vector (set in CCMModelNumExecutor)
        dx  % initialize the discritized spacing. (set in CCMModelNumExecutor)
    end
    methods
        
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
        function result = get.x(obj)
            r1 = obj.Rc;
            
            %switch to m coordinates
            m1 = (r1^3)/3;
            
            % grid spacing
            dx1 = m1/(obj.xnum);
            
            % make grid
            result = linspace(0, m1, obj.xnum);
        end
        function result = get.dx(obj)
            r1 = obj.Rc;
            
            %switch to m coordinates
            m1 = (r1^3)/3;
            
            % grid spacing
            result = m1/(obj.xnum);
            
        end
    end
end