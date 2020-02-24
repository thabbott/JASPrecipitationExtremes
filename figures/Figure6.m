%% ABC2020 Figure 6

% Setup
addpath(genpath('..\matlab'));
SAM_defineConstants;
SOM_defineConstants;

sst = [280 285 290 295 300 305];
cols = viridis(numel(sst));
sstm = 0.5*(sst(2:end) + sst(1:end-1));
get_export_fname = @(sst, fname) sprintf('../data/ch_cam%dri0/%s', ...
    sst, fname);

figure(standard_figure('half-page')); hold on;
plot(sst, 100*SAM_L_c./(SAM_R_v.*sst.^2), 'k--', 'LineWidth', 1.5);
plot(sst, zeros(size(sst)), 'k--', 'LineWidth', 1.5);
xlabel('SST (K)');
ylabel('Rate of change (% K^{-1})');

wscal = zeros(numel(sstm),2);
Mscal = zeros(numel(sstm),2);
rhoscal = zeros(numel(sstm),2);
ii = 1;
for ss = sst
    
    load(get_export_fname(ss, 'extremes9999.mat'));
    if ii < numel(sst)
        new = load(get_export_fname(sst(ii+1), 'extremes9999.mat'));
        imax = find(mass_flux_cond9999 == max(mass_flux_cond9999),1);
        new.imax = find(new.mass_flux_cond9999 == max(new.mass_flux_cond9999), 1);
        Mscal(ii,1) = mass_flux_cond9999(imax);
        rhoscal(ii,1) = rho(imax);
        wscal(ii,1) = Mscal(ii,1)/rhoscal(ii,1);
        Mscal(ii,2) = new.mass_flux_cond9999(new.imax);
        rhoscal(ii,2) = new.rho(new.imax);
        wscal(ii,2) = Mscal(ii,2)/rhoscal(ii,2);
    end
    ii = ii+1;
end

scal = @(x) 200*(x(:,2) - x(:,1))./(x(:,2) + x(:,1))./diff(sst);
plot(sstm, scal(Mscal), '-o', 'LineWidth', 1.5, 'Color', [0 0.5 0]);
plot(sstm, scal(wscal), '--o', 'LineWidth', 1.5, 'Color', [0 0.5 0]);
plot(sstm, scal(rhoscal), '-.o', 'LineWidth', 1.5, 'Color', [0 0.5 0]);
text(281, 1.4, 'Dyn. mode', 'Color', 'black');
text(284, 3.8, 'Change in w', 'Color', 'black');
text(282, -2, 'Change in \rho', 'Color', 'black');
pngandpdf('Figure6')
ylim([-3 5]);
