%% ABC2020 Figure 3

% Setup
addpath(genpath('..\matlab'));
SAM_defineConstants;
SOM_defineConstants;

sst = [280 285 290 295 300 305];
cols = viridis(numel(sst));
sstm = 0.5*(sst(2:end) + sst(1:end-1));
get_export_fname = @(sst, fname) sprintf('../data/ch_cam%dri0/%s', ...
    sst, fname);

%% 99.99th percentile
ii = 1;
nboot = 1000;
cond = zeros(numel(sst),1);
cond_boot = zeros(numel(sst), nboot);
condm = zeros(numel(sst)-1,1);
condm_ref = zeros(numel(sst)-1,1);
condm_boot = zeros(numel(sst)-1, nboot);
condm_boot_ref = zeros(numel(sst)-1, nboot);
condt = zeros(numel(sst)-1,1);
condt_ref = zeros(numel(sst)-1,1);
condt_boot = zeros(numel(sst)-1, nboot);
condt_boot_ref = zeros(numel(sst)-1, nboot);
prec = zeros(numel(sst),1);
prec_boot = zeros(numel(sst), nboot);
eff = zeros(numel(sst),1);
eff_boot = zeros(numel(sst), nboot);

for ss = sst
    
    % Compute and save column condensation
    load(get_export_fname(ss, 'extremes9999.mat'), 'mass_flux_cond9999', ...
        'precip9999', 'p', 'rho', 'dz', 'dzqn_ref');
    load(get_export_fname(ss, 'extremes9999_boot.mat'), ...
        'mass_flux_cond9999_boot');
    load(get_export_fname(ss, 'precip9999_boot.mat'));
    cond(ii) = sum(mass_flux_cond9999.*dzqn_ref.*dz);
    cond_boot(ii,:) = sum(mass_flux_cond9999_boot.*dzqn_ref.*dz, 1);

    % Save precipitation
    prec(ii) = precip9999;
    prec_boot(ii,:) = precip9999_boot;
    
    % Compute and save precipitation efficiency
    eff(ii) = prec(ii)./cond(ii);
    eff_boot(ii,:) = prec_boot(ii,:)./cond_boot(ii,:);

    % Compute condensation with dynamic and thermodynamic changes
    if ii < numel(sst)
        
        % Thermodynamic mode
        new = load(get_export_fname(sst(ii+1), 'extremes9999.mat'), ...
            'mass_flux_cond9999');
        mass_flux_mid = ...
            0.5 * (new.mass_flux_cond9999 + mass_flux_cond9999);
        new = load(get_export_fname(sst(ii+1), 'extremes9999_boot.mat'), ...
            'mass_flux_cond9999_boot');
        mass_flux_mid_boot = ...
            0.5 * (new.mass_flux_cond9999_boot + mass_flux_cond9999_boot);
        condt_ref(ii) = sum(mass_flux_mid.*dzqn_ref.*dz);
        condt_boot_ref(ii,:) = sum(mass_flux_cond9999_boot.*dzqn_ref.*dz, 1);        
        load(get_export_fname(sst(ii+1), 'extremes9999.mat'), 'dzqn_ref', ...
            'rho', 'dz');
        condt(ii) = sum(mass_flux_mid.*dzqn_ref.*dz);
        condt_boot(ii,:) = sum(mass_flux_mid_boot.*dzqn_ref.*dz, 1);
        
        load(get_export_fname(ss, 'extremes9999.mat'), 'dzqn_ref', ...
            'rho', 'dz');
        new = load(get_export_fname(sst(ii+1), 'extremes9999.mat'), 'dzqn_ref', ...
            'rho', 'dz');
        dzqn_ref_mid = 0.5 * (dzqn_ref + new.dzqn_ref);
        condm_ref(ii) = sum(mass_flux_cond9999.*dzqn_ref_mid.*dz);
        condm_boot_ref(ii,:) = sum(mass_flux_cond9999_boot.*dzqn_ref_mid.*dz, 1);
        load(get_export_fname(sst(ii+1), 'extremes9999.mat'), ...
            'mass_flux_cond9999');
        load(get_export_fname(sst(ii+1), 'extremes9999_boot.mat'), ...
            'mass_flux_cond9999_boot');
        condm(ii) = sum(mass_flux_cond9999.*dzqn_ref_mid.*dz);
        condm_boot(ii,:) = sum(mass_flux_cond9999_boot.*dzqn_ref_mid.*dz, 1);
    end
    ii = ii+1;
end

figure(standard_figure('half-page')); hold on;
plot(sst, 100*SAM_L_c./(SAM_R_v.*sst.^2), 'k--', 'LineWidth', 1.5);
plot(sst, zeros(size(sst)), 'k--', 'LineWidth', 1.5);
xlabel('SST (K)');
ylabel('Rate of change (% K^{-1})');
ylim([-5 12]);    

scal = @(x,y) 2*100*(x - y)./(x + y)./diff(sst');
sstm = 0.5*(sst(1:end-1) + sst(2:end));

condm_boot = scal(condm_boot, condm_boot_ref(1:end,:));
condm_boot_ext = [min(condm_boot, [], 2) max(condm_boot, [], 2)];
condm_boot = std(condm_boot, 0, 2);
condt_boot = scal(condt_boot, condt_boot_ref(1:end,:));
condt_boot_ext = [min(condt_boot, [], 2) max(condt_boot, [], 2)];
condt_boot = std(condt_boot, 0,  2);
cond_boot = scal(cond_boot(2:end,:), cond_boot(1:end-1,:));
cond_boot_ext = [min(cond_boot, [], 2) max(cond_boot, [], 2)];
cond_boot = std(cond_boot, 0, 2);
prec_boot = scal(prec_boot(2:end,:), prec_boot(1:end-1,:));
prec_boot_ext = [min(prec_boot, [], 2) max(prec_boot, [], 2)];
prec_boot = std(prec_boot, 0, 2);
eff_boot = scal(eff_boot(2:end,:), eff_boot(1:end-1,:));
eff_boot_ext = [min(prec_boot, [], 2) max(prec_boot, [], 2)];
eff_boot= std(eff_boot, 0, 2);


errorbar(sstm, scal(prec(2:end), prec(1:end-1)), ...
    scal(prec(2:end), prec(1:end-1)) - prec_boot_ext(:,1),...
    -(scal(prec(2:end), prec(1:end-1)) - prec_boot_ext(:,2)), ...
    '-o', ...
    'Color', 'black', 'LineWidth', 1.5);
errorbar(sstm, scal(cond(2:end), cond(1:end-1)), ...
    scal(cond(2:end), cond(1:end-1)) - cond_boot_ext(:,1),...
    -(scal(cond(2:end), cond(1:end-1)) - cond_boot_ext(:,2)), ...
    '-o', ...
    'Color', 'blue', 'LineWidth', 1.5);
errorbar(sstm, scal(eff(2:end), eff(1:end-1)), ...
    scal(eff(2:end), eff(1:end-1)) - eff_boot_ext(:,1), ...
    -(scal(eff(2:end), eff(1:end-1)) - eff_boot_ext(:,2)), ...
    '-o', ...
    'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
errorbar(sstm, scal(condt(1:end), condt_ref(1:end)), ...
    scal(condt(1:end), condt_ref(1:end)) - condt_boot_ext(:,1),...
    -(scal(condt(1:end), condt_ref(1:end)) - condt_boot_ext(:,2)), ...
    '-o', ...
    'Color', 'red', 'LineWidth', 1.5);
errorbar(sstm, scal(condm(1:end), condm_ref(1:end)), ...
    scal(condm(1:end), condm_ref(1:end)) - condm_boot_ext(:,1),...
    -(scal(condm(1:end), condm_ref(1:end)) - condm_boot_ext(:,2)), ...
    '-o', ...
    'Color', [0 0.5 0], 'LineWidth', 1.5);
text(282, 10.5, 'Condensation', 'Color', 'blue');
text(292.5, 9.8, 'Precipitation', 'Color', 'black');
text(282.5, 4.8, 'Therm. mode', 'Color', 'red');
text(285.5, 3, 'Dyn. mode', 'Color', [0 0.5 0]);
text(295, -4, 'Prec. Eff.', 'Color', [0.5 0.5 0.5]);
title('p = 99.99', 'FontWeight', 'normal');
title('(b)', 'FontWeight', 'bold');
pngandpdf('Figure3B');

%% Plot actual values and compare with CC scalings
ii = 1;
nboot = 1000;
cond = zeros(numel(sst),1);
cond_boot = zeros(numel(sst), nboot);
condm = zeros(numel(sst)-1,1);
condm_ref = zeros(numel(sst)-1,1);
condm_boot = zeros(numel(sst)-1, nboot);
condm_boot_ref = zeros(numel(sst)-1, nboot);
condt = zeros(numel(sst)-1,1);
condt_ref = zeros(numel(sst)-1,1);
condt_boot = zeros(numel(sst)-1, nboot);
condt_boot_ref = zeros(numel(sst)-1, nboot);
prec = zeros(numel(sst),1);
prec_boot = zeros(numel(sst), nboot);
eff = zeros(numel(sst),1);
eff_boot = zeros(numel(sst), nboot);
for ss = sst
    
    % Compute and save column condensation
    load(get_export_fname(ss, 'extremes9999.mat'), 'mass_flux_cond9999', ...
        'precip9999', 'p', 'rho', 'dz', 'dzqn_ref');
    load(get_export_fname(ss, 'extremes9999_boot.mat'), ...
        'mass_flux_cond9999_boot');
    load(get_export_fname(ss, 'precip9999_boot.mat'));
    cond(ii) = sum(mass_flux_cond9999.*dzqn_ref.*dz);
    cond_boot(ii,:) = sum(mass_flux_cond9999_boot.*dzqn_ref.*dz, 1);

    % Save precipitation
    prec(ii) = precip9999;
    prec_boot(ii,:) = precip9999_boot;
    
    % Compute and save precipitation efficiency
    eff(ii) = prec(ii)./cond(ii);
    eff_boot(ii,:) = prec_boot(ii,:)./cond_boot(ii,:);

    % Compute condensation with dynamic and thermodynamic changes
    if ii < numel(sst)
        
        % Thermodynamic mode
        new = load(get_export_fname(sst(ii+1), 'extremes9999.mat'), ...
            'mass_flux_cond9999');
        mass_flux_mid = ...
            0.5 * (new.mass_flux_cond9999 + mass_flux_cond9999);
        new = load(get_export_fname(sst(ii+1), 'extremes9999_boot.mat'), ...
            'mass_flux_cond9999_boot');
        mass_flux_mid_boot = ...
            0.5 * (new.mass_flux_cond9999_boot + mass_flux_cond9999_boot);
        condt_ref(ii) = sum(mass_flux_mid.*dzqn_ref.*dz);
        condt_boot_ref(ii,:) = sum(mass_flux_cond9999_boot.*dzqn_ref.*dz, 1);        
        load(get_export_fname(sst(ii+1), 'extremes9999.mat'), 'dzqn_ref', ...
            'rho', 'dz');
        condt(ii) = sum(mass_flux_mid.*dzqn_ref.*dz);
        condt_boot(ii,:) = sum(mass_flux_mid_boot.*dzqn_ref.*dz, 1);
        
        load(get_export_fname(ss, 'extremes9999.mat'), 'dzqn_ref', ...
            'rho', 'dz');
        new = load(get_export_fname(sst(ii+1), 'extremes9999.mat'), 'dzqn_ref', ...
            'rho', 'dz');
        dzqn_ref_mid = 0.5 * (dzqn_ref + new.dzqn_ref);
        condm_ref(ii) = sum(mass_flux_cond9999.*dzqn_ref_mid.*dz);
        condm_boot_ref(ii,:) = sum(mass_flux_cond9999_boot.*dzqn_ref_mid.*dz, 1);
        load(get_export_fname(sst(ii+1), 'extremes9999.mat'), ...
            'mass_flux_cond9999');
        load(get_export_fname(sst(ii+1), 'extremes9999_boot.mat'), ...
            'mass_flux_cond9999_boot');
        condm(ii) = sum(mass_flux_cond9999.*dzqn_ref_mid.*dz);
        condm_boot(ii,:) = sum(mass_flux_cond9999_boot.*dzqn_ref_mid.*dz, 1);
    end
    ii = ii+1;
end

% Define CC scaling line
ssts = linspace(sst(1), sst(end), 100);
psat = SOM_qsat(ssts, 1e5*ones(size(ssts)));
CCscal = psat/psat(1);

figure(standard_figure('half-page')); hold on;
errorbar(sst, 3600*prec, 3600*(max(prec_boot, [], 2) - prec), ...
    3600*(prec - min(prec_boot, [], 2)), ...
    '-o', 'Color', 'black', 'LineWidth', 1.5);
errorbar(sst, 3600*cond, 3600*(cond - min(cond_boot, [], 2)), ...
    3600*(max(cond_boot, [], 2) - cond), ...
    '-o', 'Color', 'blue', 'LineWidth', 1.5);
set(gca, 'YScale', 'log');
set(gca, 'YMinorTick', 'on');
rstart = 1e3;
yl = ylim;
while rstart*CCscal(end) > yl(1)
    plot(ssts, rstart*CCscal, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    rstart = rstart/2;
end
ylim(yl);
xlabel('SST (K)');
ylabel('Rate (mm hr^{-1})');
leg = legend('Precipitation', 'Condensation', 'Location', 'northwest');
leg.Box = 'off';
leg.ItemTokenSize = [10,1];
title('(a)');
pngandpdf('Figure3A');
