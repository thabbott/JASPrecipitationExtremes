%% ABC2020 Figure 2

% Setup
addpath(genpath('..\matlab'));
latex = 0;
set_plot_styling;
SAM_defineConstants;
SOM_defineConstants;

% Load SAM precipitation data
load('../data/ch_cam305ri0/exported_data.mat', 'precip');
precip = precip*3600;

% Load Nauru precipitation data
load('../data/precip_timeseries_1minavgs_nauru_org_01Apr2001_15Aug2006.mat');
nauru = precip_timeseries_1minavgs_nauru_org_01Apr2001_15Aug2006;

% Calculate PDFs
edges = [0,logspace(0,log10(max(precip(:))),round(sqrt(sum(precip(:) >= 1))))];
[spdf, sc] = pdfcounts(precip(:), edges);
edges = [0,logspace(0,log10(max(nauru(:))),round(sqrt(sum(nauru(:) >= 1))))];
[npdf, nc] = pdfcounts(nauru(:), edges);

% Load SAM percentile ranges
load('../data/ch_cam305ri0/extremes9999.mat', 'precip9999_all');
load('../data/ch_cam305ri0/extremes999.mat', 'precip999_all');

% Calculate Nauru 99.99th percentile
nauru_sort = sort(nauru(~isnan(nauru)));
ip = round(numel(nauru_sort)*0.9999);
nauru_9999 = nauru_sort(ip);

%% Plot
gray = [0.5 0.5 0.5];
figure(standard_figure('half-page')); hold on;
rectangle('Position', [3600*min(precip999_all(:)), ...
                       1e-8, ...
                       3600*(max(precip999_all(:)) - min(precip999_all(:))), ...
                       1e-3 - 1e-8], ...
          'FaceColor', gray, 'EdgeColor', 'none');
text(3600*max(precip999_all(:))/2, 2e-3, 'SAM p \approx 99.9', ...
    'FontSize', 9, 'Color', gray);
rectangle('Position', [3600*min(precip9999_all(:)), ...
               1e-8, ...
               3600*(max(precip9999_all(:)) - min(precip9999_all(:))), ...
               1e-4 - 1e-8], ...
          'FaceColor', gray, 'EdgeColor', 'none');
text(3600*max(precip9999_all(:))/2, 2e-4, 'SAM p \approx 99.99', ...
    'FontSize', 9, 'Color', gray);
plot([nauru_9999, nauru_9999], [1e-8, 1e-4 - 1e-8], ...
    'k--', 'LineWidth', 2);
text(110, 4e-5, 'Nauru p = 99.99', 'FontSize', 9, 'Color', 'black');
l = plot(sc, spdf, 'b', 'LineWidth', 1.5);
l = [l, plot(nc, npdf, 'k', 'LineWidth', 1.5)];
set(gca, 'YScale', 'log');
xlim([1, 250]);
ylim([1e-8, 1e-2])
xlabel('Precipitation rate (mm hr^{-1})');
ylabel('Probability density (hr mm^{-1}, log)');
leg = legend(l, {'SAM 305 K (instant.)', 'Nauru (1 min. avg.)'}, ...
    'Location', 'northeast');
leg.Box = 'off';
leg.ItemTokenSize = [10,1];
title('(a)');
pngandpdf('Figure2A');

%% Plot scaling of precipitation extremes at different percentiles
sst = [280, 285, 290, 295, 300, 305];
p99 = zeros(numel(sst), 1);
p999 = zeros(numel(sst), 1);
p9999 = zeros(numel(sst), 1);
p99999 = zeros(numel(sst), 1);
ii = 1;
for ss = sst
    load(sprintf('../data/ch_cam%dri0/exported_data.mat', ss), 'precip');
    precip = precip*3600;
    p = quantile(precip(:), [0.99, 0.999, 0.9999, 0.99999]);
    p99(ii) = p(1);
    p999(ii) = p(2);
    p9999(ii) = p(3);
    p99999(ii) = p(4);
    ii = ii+1;
end
cols = plasma(4);
figure(standard_figure('half-page')); hold on;
scal = @(x) 200 * ...
    (x(2:end) - x(1:end-1)) ./ (x(2:end) + x(1:end-1)) ./ diff(sst');
sstm = 0.5 * (sst(2:end) + sst(1:end-1));
plot(sst, zeros(size(sst)), 'k--', 'LineWidth', 1.5);
plot(sst, 100 * SAM_L_c ./ (SAM_R_v * sst.^2), 'k--', 'LineWidth', 1.5);
l = [...
    plot(sstm, scal(p99), '-o', 'LineWidth', 1.5, 'Color', cols(1,:)), ...
    plot(sstm, scal(p999), '-o', 'LineWidth', 1.5, 'Color', cols(2,:)), ...
    plot(sstm, scal(p9999), '-o', 'LineWidth', 1.5, 'Color', cols(3,:)), ...
    plot(sstm, scal(p99999), '-o', 'LineWidth', 1.5, 'Color', cols(4,:))];
leg = legend(l, {'p = 99', '99.9', '99.99', '99.999'}, ...
    'Position', [0.7343, 0.6723, 0.1704, 0.2565]);
leg.Box = 'off';
leg.ItemTokenSize = [10, 1];
ylim([-1, 10]);
xlabel('SST (K)');
ylabel('Rate of change (% K^{-1})');
title('(b)');
pngandpdf('Figure2B');
