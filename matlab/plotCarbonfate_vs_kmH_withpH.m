% Calculate carbon fate breakdown for varying kmH and plot bicarbonate
% transport per CO2 fixed

% Params for 3 cases.
% 1: cell has full CCM.
% 2: cell has active transport and CA, no carboxysome.
% 3: cell has just RuBisCO.
ccm_params = CCMParams_Csome;
ccm_params_cell = CCMParams_NoCsome;
ccm_params_just_rbc = CCMParams_NoCsome;

% We are sweeping over cytoplasmic pH and calculating implied
% permeabilities of the cell membrane to bicarbonate.
pH = linspace(2, 11, 30);
kmH = zeros(30);

Hmax = 30000;   % Maximum cytoplasmic bicarbonate conc. in uM
kc_opt = 1e-4;  % optimal carboxysome permeability from previous work.
alpha = 0;      % assume no conversion of cytoplasmic CO2 to bicarbonate.

% Calculate carbon fate for cells w/ carboxysomes
for i =  1:length(pH)  
    ccm_params.pH = pH(i);
    ccm_params.k = kc_opt;
    % CCMParams calculates the implied kmH for us.
    kmH(i) = ccm_params.kmH; 
    
    % Calculate the optimal jc.
    jc_opt(i) = ccm_params.CalcOptimalJc(Hmax);
    ccm_params.jc = jc_opt(i);
    
    % Run the model
    executor = FullCCMModelExecutor(ccm_params);
    res = executor.RunAnalytical();
    Hin(i) = res.Hin_pm;
    CratewO(i) = res.CratewO_pm;
    Ccyto(i) = res.c_cyto_uM;
    Ccsome(i) = res.c_csome_uM;
end

% Calculate carbon fate for cells w/o carboxysomes
for i =  1:length(pH)
    ccm_params_cell.pH = pH(i);
    
    % CCMParams calculates the implied kmH for us.
    jc_opt_cell(i) = ccm_params_cell.CalcOptimalJc(Hmax);
    ccm_params_cell.jc = jc_opt_cell(i);
    
    % note: k (ci uptake rate into cbsome) meaningless in this case as 
    % there is no no cbsome
    
    % Run the model
    executor = NoCsomeModelExecutor(ccm_params_cell);
    res = executor.RunAnalytical();
    Hin_cy(i) = res.Hin_pm;
    CratewO_cy(i) = res.CratewO_pm;
    Ccyto_cy(i) = res.c_cyto_uM;
end

% Calculate carbon fate for cells w/o just RuBisCO
% TODO(flamholz): may need to account for differt pH inside and outside
% as our approximation assumes flux in dominated by active transport.
for i =  1:length(pH)
    % set permeability kmH
    ccm_params_just_rbc.pH = pH(i);
    ccm_params_just_rbc.jc = 0;  % no active transport
    ccm_params_just_rbc.NCA = 1e-3; % no carbonic anhydrase
    % note: kc meaningless in this case as no cbsome
    % note: using > 0 NCA for now to avoid some stupid error in underlying
    % code.
    
    executor = NoCsomeModelExecutor(ccm_params_just_rbc);
    res = executor.RunAnalytical();
    Hin_just_rbc(i) = res.Hin_pm;
    CratewO_just_rbc(i) = res.CratewO_pm;
    Ccyto_just_rbc(i) = res.c_cyto_uM;
end

% calculate pH from the kmH's selected
cost_h = Hin./CratewO;
cost_cy_h = Hin_cy./CratewO_cy;
cost_just_rbc_h = Hin_just_rbc./CratewO_just_rbc;

figure(1)
loglog(kmH, cost_h, 'Color', 'k', 'LineWidth', 3);
hold on
text(10*kmH(end), 1.7*cost_h(end), 'CCM cost', 'Color', 'r');

loglog(kmH, cost_cy_h, 'Color', 'k', 'LineWidth', 3);
plot([1e-7 1e-7], [1e-1 2e7], 'Color', 'm', 'LineWidth', 20)
text(1.7e-7, 5e4, 'HCO_3^- permeability', 'Color', 'm')

plot([3e-3 3e-3], [1e-1 2e7], 'Color', 'c', 'LineWidth', 20)
text(7e-5, 5e4, 'H_2CO_3 permeability', 'Color', 'c')

% 3 ATP, 2 NADPH per CO2 fixed in Calvin Cycle.
% 4 H+ per ATP, 10 H+ per NADPH. 
cost_fixation = 4*3.0 + 10*2.0;
plot([kmH(1) kmH(end)], [cost_fixation cost_fixation], '--b');
%text(kmH(1)/1100.0, 50, 'Calvin Cycle Fixation', 'Color', 'b');

plot([3e-4 3e-4], [1e-1 3e3], '--k');
plot([kmH(1) kmH(end)], [3e3 3e3], '--k');
%text(kmH(1)/1100.0, 5e3, 'Original Model', 'Color', 'k');

% cost is 1 nadph, 1 atp and 1 co2 times specificity ratio
% when RuBisCO is saturated: S = vc/vo ~ 13;
% not saturated: S = (vc/Kc)/(vo/Ko)~ 50 %SavirMilo paper
% Is rubisco saturated?
frac_oxygenation_cell = ccm_params.O ./ (Ccyto_cy .* ccm_params.S_sat);
cost_2pg_cell_h = frac_oxygenation_cell .* (10 + 4 + 1.5*4 + 10);
total_cost_cell_h = cost_2pg_cell_h + cost_cy_h + cost_fixation;

plot(kmH, total_cost_cell_h, '--g');
text(kmH(1)/1100.0, 1.5*total_cost_cell_h(end), 'No CCM', 'Color', 'g');

frac_oxygenation_noca = ccm_params.O ./ (Ccyto_just_rbc .* ccm_params.S_sat);
cost_2pg_noca_h = frac_oxygenation_noca .* (10 + 4 + 1.5*4 + 10);
total_cost_noca_h = cost_2pg_noca_h + cost_just_rbc_h + cost_fixation;

plot(kmH, total_cost_noca_h, '--m');
%text(kmH(1)/1100.0, 1.5*total_cost_noca_h(end), 'Just RuBisCO', 'Color', 'g');

haxes1 = gca; % handle to axes
set(haxes1,'XColor','k','YColor','k')
axis([kmH(end) kmH(1) 1e-1 2e7])
xlabel('k_m^H, Effective Membrane Permeability of HCO_3^* (cm/s)')
haxes1_pos = get(haxes1,'Position'); % store position of first axes
haxes2 = axes('Position',haxes1_pos,...
              'XAxisLocation','top',...
              'YAxisLocation','left',...
              'Yscale', 'log', ...
              'Color','none');
hold on
loglog(pH, cost_h,'Parent', haxes2,'Color','r', 'LineWidth', 3);
set(haxes2, 'XDir','Reverse')
axis([pH(1) pH(end) 1e-1 2e7])
ylabel('Energetic Cost H^+ / (CO_2 fixed)')
xlabel('Cytoplasmic pH')

