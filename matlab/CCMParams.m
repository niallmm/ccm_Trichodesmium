classdef CCMParams
    % Object defining CCM parameters - encapsulates various dependent
    % calculations of rates and volumes.
    
    % Mutable properties that may depend on the model context.
    properties
        jc = 0.6;        % active uptake rate of HCO3- (cm/s)
        kcC = 1e-4;      % permeability of carboxysome to CO2 (cm/s)
        kcH = 1e-4;      % permeability of carboxysome to HCO3- cm/s
        Rc = 5e-6;       % radius of c-some (cm)
        Rb = 5e-5;       % radius of cell (cm)
        D = 1e-5;        % diffusion constant (cm^2/s)
        
        kmC = 0.3;       % cm/s permiability of outer membrane to CO2
        alpha = 0;       % reaction rate of conversion of CO2 to HCO3- at the cell membrane (cm/s)
        Cout = 15;       % microM from Henry's law.
        Hout = 1900; 
        
        % values at pH 7.8 -- will be used to re-scale pH dependence in
        % function
        % kRub, Km_8, KO, Ssat are all for Synechococcus 6301 from Savir
        kRub = 11.6;     % rxns/s maximum reaction rate at single active site
        NRub = 2000;     % number of RuBisCO active sites
        pHmeas = 7.8;    % pH at which the RuBisCO parameters were measured
        Km_meas = 340;   % half max reaction rate of RuBisCO, uM
        KO_meas = 972;   % uM
        O = 300;         % uM, calculated from ambient O2
        S_sat = 43;      % Specificity ratio when RuBisCO is saturated
        
        kCAH = 4.6e4;    % maximum rate of bicarbonate dehydration by CA /active site @ pH = 8
        kCAC = 8e4;      % maximum rate of carbon dioxide hydration by CA /active site
        NCA = 100;       % number of active sites
        Kca = 3.2*1e3;   % uM half max reaction rate for carbon dioxide hydration
        Kba = 9.3*1e3;   % uM half max reaction rate for bicarbonate dehydration
        
        pH = 8;          % pH in cytosol
        pH_csome = 8;    % pH in csome
        kmH_base = 3e-3;      % cm/s this is the permeability of the membrane to pure H2CO3
        h_cyto_exp = 30000;   % uM of inorganic carbon expected in the cytosol
        
        pH_out = 7;          % extracellular pH, 7 is around freshwater
        I_out = 0.05;        % ionic strength of solution:
                             % 0.2 is representative of cytosol
        I_in =  0.2;         % 0.05 is representative of freshwater

        salt = 0;           % if freshwater (salt =0) we calculate pKas analytically
                            % if saltwater (salt = 1) we set pKas empirically
        pHoff= 0;           % flag to turn the pH dependence off when running numerical simulations.
        
    end
    
    % Values that cannot be edited by client code since they are physical
    % constants.
    properties (Constant)
        Na = 6.022e23;       % Avogadro's number is constant, of course.
        RT = 2.4788;         % (R = 8.314e-3 kJ/(K*mol))*(298.15 K)
        
        % Free energies of formation, charges and proton counts. 
        CO2_delG0 = -623.2;
        CO2_NH = 2;
        CO2_z = 0;
        H2CO3_delG0 = -606.3;
        H2CO3_NH = 2;
        H2CO3_z = 0;
        HCO3_delG0 = -586.8;
        HCO3_NH = 1;
        HCO3_z = -1;
        CO3_delG0 = -527.8;
        CO3_NH = 0;
        CO3_z = -2;
    end
    
    properties (Dependent)
        % volumes and pH dependent RuBisCO reaction rates
        Vcell   % volume of cell
        Vcsome  % volume of carboxysome
        SAcell  % surface area of the cell
        
        % Dependent paramters for the case that CA & RuBisCO are uniformly
        % co-localized to the carboxysome.
        VmaxCsome    % uM/s RuBisCO max reaction rate/concentration
        
        % needed to calculate the oxygen reaction rate
        VbaCsome     % maximum rate of bicarbonate dehydration by CA
        VcaCsome     % maximum rate of carbon dioxide hydration by CA
        
        % Dependent paramters for the case that CA & RuBisCO are uniformly
        % distributed through the cytoplasm.
        VmaxCell     % uM/s RuBisCO max reaction rate/concentration
        VbaCell      % maximum rate of bicarbonate dehydration by CA
        VcaCell      % maximum rate of carbon dioxide hydration by CA
        
        kRub_pH     % pH dependent RuBisCO reaction rate/s at single reaction site
        Km          % pH dependent RuBisCO 1/2 max concentration
        KO          % pH dependent RuBisCO 1/2 max concentration for oxygen,
                    % interpreted from pH dependence of Km
    end
    methods
        % calculating reaction volumes and pH dependent RuBisCOreaction rates
        function value = get.Vcell(obj)
            value = 4*pi*obj.Rb^3/3;
        end
        function value = get.Vcsome(obj)
            value = 4*pi*obj.Rc^3/3;
        end
        
        function value = get.SAcell(obj)
            value = 4*pi*obj.Rb^2;
        end
        function value = get.kRub_pH(obj)
            if obj.pHoff ==1
                if obj.pH_csome ~=8
                    obj.pHoff
                    warning('kRub: pH dependence is turned off but you requested a pH other than 8')
                end
                value = 10.4231;
            else
                rxn = load('pH_Vmax.mat');
                VmaxpHmeas = interp1(rxn.pH, rxn.Vmax, obj.pHmeas, 'pchip'); % find value of Vmax at pH 8 from data
                if (obj.pH_csome > 8.36) || (obj.pH_csome < 6)
                    warning('pH is out of range of experimental RuBisCO measurements (6 to 8.46) but we are extrapolating')
                    warning('we are extrapolating out to pH 8.9')
                end
                if (obj.pH_csome > 8.9) || (obj.pH_csome < 6)
                    
                    if (obj.pH_csome < 6)
                        temp = interp1(rxn.pH, rxn.Vmax, 6, 'pchip'); % find value of Vmax at pH we want from data
                        value = temp*obj.kRub/VmaxpHmeas;
                    end
                    if (obj.pH_csome>8.9)
                        temp = interp1(rxn.pH, rxn.Vmax, 8.9, 'pchip');
                        value = temp*obj.kRub/VmaxpHmeas;
                    end
                else
                    temp = interp1(rxn.pH, rxn.Vmax, obj.pH_csome, 'pchip'); % find value of Vmax at pH we want from data
                    value = temp*obj.kRub/VmaxpHmeas; % scale the interpolated value by our known kRub rate at pH 8
                    % need to scale because the data was not in the right units.
                    % now it is scaled to rxns/s per active site
                end
            end
        end
        function value = get.Km(obj)
            if obj.pHoff ==1
                if obj.pH_csome~=8
                    obj.pHoff
                    warning('Km: pH dependence is turned off but you requested a pH other than 8')
                end
                value = 276.9778;
            else
                rxn = load('pH_Km.mat');
                KmpHmeas = interp1(rxn.pH, rxn.Km, obj.pHmeas, 'pchip'); % find value of Vmax at pH 8 from data
                if (obj.pH_csome > 8.46) || (obj.pH_csome < 6)
                    warning('pH is out of range of experimental RuBisCO measurements (6 to 8.46) but we are extrapolating')
                    warning('we can extrapolate out to pH 8.9')
                end
                if (obj.pH_csome > 8.9) || (obj.pH_csome < 6)
                    
                    if (obj.pH_csome <6)
                        value = interp1(rxn.pH, rxn.Km*obj.Km_meas/KmpHmeas, 6, 'pchip'); % Km was in mM
                    end
                    if obj.pH_csome >8.9
                        value = interp1(rxn.pH, rxn.Km*obj.Km_meas/KmpHmeas, 8.9, 'pchip');
                    end
                else
                    value = interp1(rxn.pH, rxn.Km*obj.Km_meas/KmpHmeas, obj.pH_csome,'pchip'); % Km was in mM
                end
            end
        end
        
        function value = get.KO(obj)
            if obj.pHoff ==1
                if obj.pH_csome ~=8
                    obj.pHoff
                    warning('KO:pH dependence is turned off but you requested a pH other than 8')
                end
                value = 791.8306;
            else
                rxn = load('pH_Km.mat');
                KmpHmeas = interp1(rxn.pH, rxn.Km, obj.pHmeas, 'pchip'); % find value of Vmax at pH 8 from data
                if (obj.pH_csome > 8.5) || (obj.pH_csome < 6)
                    warning('pH is out of range of experimental RuBisCO measurements (6 to 8.46) but we are extrapolating')
                    warning('we are extrapolating out to pH 8.9')
                end
                if (obj.pH_csome > 8.9) || (obj.pH_csome < 6)
                    
                    if (obj.pH_csome <6)
                        value = interp1(rxn.pH, rxn.Km*obj.KO_meas/KmpHmeas, 6, 'pchip'); % Km was in mM
                    end
                    if obj.pH_csome>8.9
                        value = interp1(rxn.pH, rxn.Km*obj.KO_meas/KmpHmeas, 8.9, 'pchip');
                    end
                else
                    value = interp1(rxn.pH, rxn.Km*obj.KO_meas/KmpHmeas, obj.pH_csome, 'pchip'); % Km was in mM
                end
            end
        end
        
        function value = get.VmaxCsome(obj)
            value = obj.kRub_pH * obj.NRub*1e6/(obj.Vcsome * obj.Na * 1e-3);
        end
        function value = get.VbaCsome(obj)
            value = obj.kCAH * obj.NCA * 1e6/(obj.Vcsome * obj.Na * 1e-3);
        end
        function value = get.VcaCsome(obj)
            value = obj.kCAC * obj.NCA * 1e6/(obj.Vcsome * obj.Na * 1e-3);
        end
        
        function value = get.VmaxCell(obj)
            value = obj.kRub_pH * obj.NRub*1e6/(obj.Vcell * obj.Na * 1e-3);
        end
        function value = get.VbaCell(obj)
            value = obj.kCAH * obj.NCA * 1e6/(obj.Vcell * obj.Na * 1e-3);
        end
        function value = get.VcaCell(obj)
            value = obj.kCAC * obj.NCA * 1e6/(obj.Vcell * obj.Na * 1e-3);
        end
    end
    properties (Dependent)
        % Gibbs free energies corrected for ionic strength of solution
        delG0I_CO2_out     % CO2 ionic corrected Gibbs free energy of formation
        delG0I_H2CO3_out   % H2CO3 ionic corrected Gibbs free energy of formation
        delG0I_HCO3_out   % HCO3- ionic corrected Gibbs free energy of formation
        delG0I_CO3_out   % CO3(2-) ionic corrected Gibbs free energy of formation
        
        delG0I_CO2_in     % CO2 ionic corrected Gibbs free energy of formation
        delG0I_H2CO3_in    % H2CO3 ionic corrected Gibbs free energy of formation
        delG0I_HCO3_in   % HCO3- ionic corrected Gibbs free energy of formation
        delG0I_CO3_in   % CO3(2-) ionic corrected Gibbs free energy of formation
        
        % pKa's for Ci pool in the cytosol, calculated assuming I = 0.2
        pKa1_cyto       % pKa for the H2CO3 to HCO3- equilibrium
        pKa2_cyto       % pKa for CO3(-2) to HCO3- equilibrium
        pKa_eff_cyto    % effective pKa for CO2 to HCO3- equilibrium
        
        % pKa's for Ci pool ouside the cell, can be for either freshwater or
        % saltwater
        pKa1_out
        pKa2_out
        pKa_eff_out
        
%         Hout   % external HCO3- concentration, calculated from
        % Cout < ---- now in first properties section (10 uM from Henry's
        % law)
        
        %pH dependent things
        kmH_in         % cm/s permiability of outer membrane to HCO3- ph dependent
        kmH_out     % same as above, but with external
        Keq         % pH dependent equilibrium constant of carbonic anhydrase
        
        
    end
    methods
        function value = get.delG0I_CO2_out(obj)
            value = obj.CO2_delG0 - 2.91482*(obj.CO2_z^2 -obj.CO2_NH)*sqrt(obj.I_out)./(1+1.6*sqrt(obj.I_out));
        end
        function value = get.delG0I_H2CO3_out(obj)
            value = obj.H2CO3_delG0 - 2.91482*(obj.H2CO3_z^2 -obj.H2CO3_NH)*sqrt(obj.I_out)./(1+1.6*sqrt(obj.I_out));
        end
        function value = get.delG0I_HCO3_out(obj)
            value = obj.HCO3_delG0 - 2.91482*(obj.HCO3_z^2 -obj.HCO3_NH)*sqrt(obj.I_out)./(1+1.6*sqrt(obj.I_out));
        end
        function value = get.delG0I_CO3_out(obj)
            value = obj.CO3_delG0 - 2.91482*(obj.CO3_z^2 -obj.CO3_NH)*sqrt(obj.I_out)./(1+1.6*sqrt(obj.I_out));
        end
        function value = get.delG0I_CO2_in(obj)
            value = obj.CO2_delG0 - 2.91482*(obj.CO2_z^2 -obj.CO2_NH)*sqrt(obj.I_in)./(1+1.6*sqrt(obj.I_in));
        end
        function value = get.delG0I_H2CO3_in(obj)
            value = obj.H2CO3_delG0 - 2.91482*(obj.H2CO3_z^2 -obj.H2CO3_NH)*sqrt(obj.I_in)./(1+1.6*sqrt(obj.I_in));
        end
        function value = get.delG0I_HCO3_in(obj)
            value = obj.HCO3_delG0 - 2.91482*(obj.HCO3_z^2 -obj.HCO3_NH)*sqrt(obj.I_in)./(1+1.6*sqrt(obj.I_in));
        end
        function value = get.delG0I_CO3_in(obj)
            value = obj.CO3_delG0 - 2.91482*(obj.CO3_z^2 -obj.CO3_NH)*sqrt(obj.I_in)./(1+1.6*sqrt(obj.I_in));
        end
        function value = get.pKa1_out(obj)
            if obj.salt == 0
                value = (obj.delG0I_HCO3_out - obj.delG0I_H2CO3_out)/(obj.RT*log(10));
            elseif obj.salt == 1
                value = 3.1;
            else
                warning('salt is set to a value other than 0 or 1 in CCMParams.m')
            end
        end
        function value = get.pKa2_out(obj)
            if obj.salt == 0
                value = (obj.delG0I_CO3_out - obj.delG0I_HCO3_out)/(obj.RT*log(10));
            elseif obj.salt == 1
                value = 8.91;
            else
                warning('salt is set to a value other than 0 or 1 in CCMParams.m')
            end
        end
        function value = get.pKa_eff_out(obj)
            if obj.salt == 0
                value = (obj.delG0I_HCO3_out - obj.delG0I_CO2_out)/(obj.RT*log(10));
            elseif obj.salt == 1
%                  value = 5.86;
                    value = 5.84;
            else
                warning('salt is set to a value other than 0 or 1 in CCMParams.m')
            end
        end
        function value = get.pKa1_cyto(obj)
            value = (obj.delG0I_HCO3_in - obj.delG0I_H2CO3_in)/(obj.RT*log(10));
        end
        function value = get.pKa2_cyto(obj)
            value = (obj.delG0I_CO3_in - obj.delG0I_HCO3_in)/(obj.RT*log(10));
        end
        function value = get.pKa_eff_cyto(obj)
            value = (obj.delG0I_HCO3_in - obj.delG0I_CO2_in)/(obj.RT*log(10));
        end
        
        
        function value = get.kmH_in(obj)
            value = obj.kmH_base*10^(obj.pKa1_cyto - obj.pH);
        end
        function value = get.kmH_out(obj)
            value = obj.kmH_base*10^(obj.pKa1_out - obj.pH_out);
        end
        
        function value = get.Keq(obj)
            if obj.pHoff == 1
                value = obj.Vca*obj.Kba/(obj.Vba*obj.Kca);
            else
                value = 10^(-obj.pKa_eff_cyto +obj.pH_csome);
            end
        end
        
%         function value = get.Hout(obj)
% %             Citotal = obj.Cout*(1+ 10^(-obj.pKa_eff_out)*(10^(obj.pKa1_out)...
% %                 +10^(-obj.pKa2_out + 2*obj.pH_out) + 10^(+obj.pH_out)));
% %             value = Citotal/(1+10^(obj.pKa_eff_out-obj.pH_out) ...
% %                 + 10^(-obj.pKa2_out+obj.pH_out) + 10^(obj.pKa1_out-obj.pH_out));
%             value = obj.Cout*10^(-obj.pKa_eff_out+obj.pH_out);
%         end
        
    end
        
        properties (Abstract)         % Non-dimensional params
            
            % See supplementary material and pdf document NonDimEqns2 for how these
            % Formulas contain the cytosol solutions.
            xi      % ratio of rate of diffusion across cell to rate of
            %dehydration reaction of carbonic anhydrase (D/Rb^2)/(Vba/Kba)
            gamma   % ratio of forward and backward carbonic anhydrase rates (Vba/Vca)
            kappa   % ratio of 1/2 max forward and 1/2 max backward carbonic anhydrase rates
            beta_c  % dc/d\rho = beta_c*c - epsilon_c
            beta_c2 % dh/d\rho = beta_h*h - epsilon_h - beta_c2*c
            beta_h  % dh/d\rho = beta_h*h - epsilon_h - beta_c2*c
            epsilon_c % dc/dp = beta_c*c - eps_c
            epsilon_h % dh/d\rho = beta_h*h - epsilon_h - beta_c2*c
            GC  % grouped params = D/(Rc^2 kcC) + 1/Rc - 1/Rb [1/cm]
            GH  % grouped params = D/(Rc^2 kcH) + 1/Rc - 1/Rb [1/cm]
            
            % Calculated appropriate to the volume in which the enzymes are
            % contained which depends on the situation (in cbsome or not).
            Vmax    % uM/s RuBisCO max reaction rate/concentration
            Vba     % maximum rate of bicarbonate dehydration by CA
            Vca     % maximum rate of carbon dioxide hydration by CA
        end
        
        methods (Abstract)
            % Calculate the optimal jc value (active bicarbonate influx).
            % Calculated differently depending on the model.
            CalcOptimalJc(obj)
        end
    end
