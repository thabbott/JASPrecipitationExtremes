%% ABC2020 reproduce panel in Figure 5
%% with alternate simulations
%% for Figure A2

% Setup
addpath(genpath('..\matlab'));
latex = 1;
set_plot_styling;
SAM_defineConstants;
SOM_defineConstants;

sst = [280 285 290 295 300 305];
cols = viridis(numel(sst));
sstm = 0.5*(sst(2:end) + sst(1:end-1));
get_export_fname = @(sst, fname) sprintf('../data/RCEMIP_SST%d_M2005_TKE15/%s', ...
    sst, fname);

%% 99.99th percentile

% Make a set of plots for each mass flux set
f1 = standard_figure('quarter-page');
f2 = standard_figure('quarter-page');
f3 = standard_figure('quarter-page');
f4 = standard_figure('quarter-page');
f5 = standard_figure('quarter-page');
f6 = standard_figure('quarter-page');
ii = 1;
condt = zeros(numel(sst)-1,1);
condb = zeros(numel(sst)-1,1);
l = [];
lbl = {};
for ss = sst

    load(get_export_fname(ss, 'extremes9999.mat'));
    load(get_export_fname(ss, 'exported_convective_w_and_buoy.mat'), ...
        'w_9999', 'rho');
    qcum = cumsum(dzqn_ref.*dz);
    qtilde = 1 - qcum/qcum(end);
    figure(f1); hold on;
    plot(mass_flux_cond9999/max(mass_flux_cond9999), ...
        p, 'LineWidth', 1.5, 'Color', cols(ii,:));
    figure(f2); hold on;
    plot(mass_flux_cond9999/max(mass_flux_cond9999), ...
        qtilde, 'LineWidth', 1.5, 'Color', cols(ii,:));
    figure(f3); hold on;
    plot(mass_flux_homo9999/max(mass_flux_homo9999), ...
        p, 'LineWidth', 1.5, 'Color', cols(ii,:));
    figure(f4); hold on;
    plot(mass_flux_homo9999/max(mass_flux_homo9999), ...
        qtilde, 'LineWidth', 1.5, 'Color', cols(ii,:));
    figure(f5); hold on;
    plot(-SAM_g * rho' .* w_9999, ...
        p, 'LineWidth', 1.5, 'Color', cols(ii,:));
    figure(f6); hold on;
    plot(-SAM_g * rho' .* w_9999, ...
        qtilde, 'LineWidth', 1.5, 'Color', cols(ii,:));
    lbl{ii} = sprintf('%dK', ss);
    ii = ii+1;
    
end
figure(f1)
set(gca, 'YDir', 'reverse');
xlim([0 1]);
ylim([0 1000]);
ylabel('Pressure (hPa)');
xlabel('$\omega$ (norm.)');
leg = legend(gca, lbl, 'Units', 'normalized', 'Position', [0.8, 0.74, 0.1, 0.1]);
leg.ItemTokenSize = [10,1];
leg.Box = 'off';
title('(a)');
pngandpdf('Figure5_a');
figure(f2);
set(gca, 'YDir', 'reverse');
xlim([0 1]);
ylim([0 1]);
xlabel('$\omega$ (norm.)');
ylabel('$\tilde{q}$ (nondim.)');
title('(d) M2005 + 1.5-TKE');
pngandpdf('Figure5_b');
figure(f3)
set(gca, 'YDir', 'reverse');
xlim([0 1]);
ylim([0 1000]);
ylabel('Pressure (hPa)');
xlabel('$\omega$ (norm.)');
title('(c)');
pngandpdf('Figure5_c');
figure(f4)
set(gca, 'YDir', 'reverse');
xlim([0 1]);
ylim([0 1]);
xlabel('$\omega$ (norm.)');
ylabel('$\tilde{q}$ (nondim.)');
title('(d)');
pngandpdf('Figure5_d');
figure(f5);
set(gca, 'YDir', 'reverse');
ylim([0 1000]);
ylabel('Pressure (hPa)');
xlabel('$\omega_{99.99}$ (norm.)');
pngandpdf('Figure5_e');
figure(f6);
set(gca, 'YDir', 'reverse');
ylim([0 1]);
xlabel('$\omega_{99.99}$ (norm.)');
ylabel('$\tilde{q}$ (nondim.)');
pngandpdf('Figure5_f');
