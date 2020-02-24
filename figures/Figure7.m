%% ABC2020 Figure 7

% Setup
addpath(genpath('..\matlab'));
latex = 1;
set_plot_styling;
SAM_defineConstants;
SOM_defineConstants;
sst = [280 285 290 295 300 305];
get_export_fname = @(sst, fname) sprintf('../Data/ch_cam%dri0/%s', ...
    sst, fname);

%% Create profiles
ps = 100000; % Surface pressure, Pa
p = linspace(ps, 1, 100); % Pressure profile, Pa

%% Compute scaling (w/ liquid CC thermodynamics)
SSTs = linspace(280, 305, 50);
dTs = 2;
qv0 = 0.81*SAM_qsatWater(SSTs-dTs, ps); % Using 0.81 b/c contouring is weird with 0.8...
wrangle = @(x) repmat(reshape(x, [1 1 numel(x)]), [size(SSTs) 1]);
[T,qv,dpqv,qstar] = CC_pseudoAdiabat(SSTs-dTs, qv0, wrangle(p), 1);
T = squeeze(T); qv = squeeze(qv); dpqv = squeeze(dpqv); qstar = squeeze(qstar);

%% Re-grid density
rho = p./(287.04.*T);
qtilde = qstar./qstar(:,1);
qtgrid = linspace(0,1,100);
rhogrid = zeros(numel(qtgrid), numel(SSTs));
for ii = 1:numel(SSTs)
    rhogrid(:,ii) = interp1(qtilde(ii,:), rho(ii,:), qtgrid);
end
    
%% Plot density
figure(standard_figure('half-page')); hold on;
contour(SSTs, qtgrid, rhogrid, 0.0:0.02:1.6, 'Color', [0.5 0.5 0.5]);
[C, h] = contour(SSTs, qtgrid, rhogrid, 0.0:0.1:1.6, 'Color', 'black');
for ss = sst
    load(get_export_fname(ss, 'extremes9999.mat'));
    qtilde_sim = 1 - cumsum(dzqn_ref.*dz)/sum(dzqn_ref.*dz);
    imax = find(mass_flux_cond9999 == max(mass_flux_cond9999),1);
    plot(ss, qtilde_sim(imax), '*', 'Color', 'black', 'LineWidth', 2);
end
title('(a) $\rho$ (kg m$^{-3}$)', 'FontWeight', 'normal');
clabel(C,h,'manual');
set(gca, 'YDir', 'reverse');
xlabel('SST (K)');
ylabel('$\tilde{q}$ (nondim.)');
pngandpdf('Figure7A');

%% Plot changes in density
dlnrhodT = 200*(rhogrid(:,2:end) - rhogrid(:,1:end-1))./...
    (rhogrid(:,2:end) + rhogrid(:,1:end-1))./diff(SSTs);
SSTm = linspace(280, 305, numel(SSTs)-1);
figure(standard_figure('half-page')); hold on;
contour(SSTm, qtgrid, dlnrhodT, -6:0.1:0, 'Color', [0.5 0.5 0.5]);
[C,h] = contour(SSTm, qtgrid, dlnrhodT, -6:0.5:0, 'Color', 'black');
title('(b) ${\rm d}\ln(\rho)/{\rm d}T$ (\% K$^{-1}$)', 'FontWeight', 'normal');
clabel(C,h, 'manual');
set(gca, 'YDir', 'reverse');
xlabel('SST (K)');
ylabel('$\tilde{q}$ (nondim.)');
pngandpdf('Figure7B');
