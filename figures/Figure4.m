%% ABC2020 Figure 4

% Setup
addpath(genpath('..\matlab'));
SAM_defineConstants;
SOM_defineConstants;

sst = [280 285 290 295 300 305];
cols = viridis(numel(sst));
alpha = 0.5;
sstm = 0.5*(sst(2:end) + sst(1:end-1));
get_export_fname = @(sst, fname) sprintf('../data/ch_cam%dri0/%s', ...
    sst, fname);

f1 = standard_figure('quarter-page');
f2 = standard_figure('quarter-page');
f3 = standard_figure('half-page');

figure(f3); hold on;
plot(sst, 100*SAM_L_c./(SAM_R_v.*sst.^2), 'k--', 'LineWidth', 1.5);
plot(sst, 140*SAM_L_c./(SAM_R_v.*sst.^2), '--', ...
    'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);
plot(sst, 65*SAM_L_c./(SAM_R_v.*sst.^2), '--', ...
    'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);
xlabel('SST (K)');
ylabel('Therm. mode (% K^{-1})');

cond1 = zeros(numel(sst)-1,1);
cond2 = zeros(numel(sst)-1,1);
lbl = {};
lbl2 = {};
l = [];
ii = 1;
for ss = sst
    
    load(get_export_fname(ss, 'extremes9999.mat'), ...
        'mass_flux_cond9999', 'p', 'dz', 'dzqn_ref', 'rho');
    
    figure(f1); hold on;
    plot(SAM_g*mass_flux_cond9999, p, ...
        'LineWidth', 1.5, 'Color', cols(ii,:));
    
    figure(f2); hold on;
    iend = find(dzqn_ref == min(dzqn_ref), 1);
    plot(1e5*dzqn_ref ./ (rho * SAM_g), p, ...
        'LineWidth', 1.5, 'Color', cols(ii,:));
    
    jj = 1;
    cond = zeros(numel(sst),1);
    for ss2 = sst        
        load(get_export_fname(ss2, 'extremes9999.mat'), ...
            'dzqn_ref', 'dz');       
        cond(jj) = sum(mass_flux_cond9999.*dzqn_ref.*dz);
        jj = jj+1;
    end
        
    figure(f3); hold on;
    l = [l, plot(sstm, ...
        200*(cond(2:end) - cond(1:end-1))./(cond(2:end) + cond(1:end-1)) ...
        ./diff(sst'), '-o', 'LineWidth', 1, 'Color', [cols(ii,:), alpha])];
    
    if ii < numel(sst)
        load(get_export_fname(sst(ii+1), 'extremes9999.mat'), ...
            'mass_flux_cond9999', 'dzqn_ref', 'dz');
        m2 = load(get_export_fname(sst(ii), 'extremes9999.mat'), ...
            'mass_flux_cond9999', 'dzqn_ref', 'dz');
        mf_mid = 0.5 * (mass_flux_cond9999 + m2.mass_flux_cond9999);
        cond1(ii) = sum(mf_mid.*dzqn_ref.*dz);
        cond2(ii) = sum(mf_mid.*m2.dzqn_ref.*m2.dz);
    end
    lbl2{ii} = sprintf('%dK', ss);
    lbl{ii+1} = [lbl2{ii}, ' K \omega'];
    ii = ii+1;
end
lbl{1} = 'True therm. mode';

figure(f3); hold on;
l = [plot(sstm, 200*(cond1 - cond2)./(cond1 + cond2)./diff(sst'), '-o', ...
    'LineWidth', 1.5, 'Color', 'red'), l];
ylim([-5 12]);
leg = gridLegend(l, 2, lbl, 'Units', 'normalized', ...
    'Location', 'southeast');
leg.Box = 'off';
leg.ItemTokenSize = [10,1];
disp('Position legend manually');
pause;
title('(c)');
pngandpdf('Figure4C');

figure(f1);
ylabel('Pressure (hPa)');
set(gca, 'YDir', 'reverse');
ylim([0 1000]);
xlim([0 65]);
xlabel('-\omega (Pa s^{-1})');
leg = ...
    legend(lbl2, 'Units', 'normalized', 'Position', [0.8, 0.73, 0.1, 0.1]);
leg.Box = 'off';
leg.ItemTokenSize = [10,1];
title('(a)');
pngandpdf('Figure4A');

figure(f2);
set(gca, 'YDir', 'reverse');
ylim([0 1000]);
xlim([0 0.03]);
xlabel('dq*/dp (g/kg/hPa)');
title('(b)');
pngandpdf('Figure4B');
