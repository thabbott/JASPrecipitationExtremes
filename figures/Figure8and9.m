%% ABC2020 Figures 8 and 9, plus a panel for Figure A4

%% Setup
addpath(genpath('../matlab'));
latex = 0;
set_plot_styling;
SAM_defineConstants;
SOM_defineConstants;

% Calculate profiles
sst = 280:5:305;        % SST
RH = 0.9;              % Free-tropospheric RH
SRH = 0.9;              % Surface RH
fnamet = ['../Data/ch_cam%dri0/%s'];
zgrid = 50:50:20e3;          % Vertical grid (m)
eps = 0.15e-3; % Entrainment rate (1/m)
fname = sprintf('LBmatch_RH_%.2f_eps_%.2f', RH, eps * 1e3);
fname = [fname, '_%s'];

f1 = standard_figure('half-page');
f2 = standard_figure('half-page');
f3 = standard_figure('half-page');
f4 = standard_figure('half-page');
f5 = standard_figure('half-page');
f6 = standard_figure('half-page');
l1 = [];
l2 = [];
l3 = [];
l4 = [];
l5 = [];
cols = viridis(numel(sst));
ii = 1;
lbl = {};

rho_qtilde = zeros(numel(sst),1);
buoyint_model_qtilde = zeros(numel(sst),1);
buoyint_sim_qtilde = zeros(numel(sst),1);
buoyint_model_220K = zeros(numel(sst),1);
buoyint_sim_220K = zeros(numel(sst),1);
w_reda_qtilde = zeros(numel(sst), 1);
w_sticky_qtilde = zeros(numel(sst), 1);
w_full_qtilde = zeros(numel(sst), 1);
w_cond_qtilde = zeros(numel(sst), 1);
w_9999_qtilde = zeros(numel(sst), 1);
w_9999_max = zeros(numel(sst), 1);
w_cond_peak = zeros(numel(sst), 1);
Tzb_all = zeros(numel(sst), 1);
D_p = zeros(numel(sst), 1);
D_q = zeros(numel(sst), 1);
D_z = zeros(numel(sst), 1);
buoy_anom = zeros(numel(sst), 1);
for ss = sst

    % Load profiles
    load(sprintf(fnamet, ss, 'extremes9999.mat'), ...
        'mass_flux_cond9999', 'tabs_ref', 'rho', 'z', 'p', ...
        'mass_flux_cond9999_all', 'dzqn_ref', 'dz');
    load(sprintf(fnamet, ss, 'exported_data.mat'), ...
        'qv_ref');
    
    % Load 99.99th percentile buoyancy from simulations
    masked = load(sprintf(fnamet, ss, ...
        'exported_convective_w_and_buoy.mat'), 'w_9999', 'buoy_9999', 'z');

    % Interpolate onto high-resolution vertical grid
    tabs_ref = interp1(z, tabs_ref, zgrid);
    mass_flux_cond9999 = interp1(z, mass_flux_cond9999, zgrid);
    mass_flux_cond9999_all = interp1(z, mass_flux_cond9999_all, zgrid);
    masked_buoy9999 = interp1(masked.z, masked.buoy_9999, zgrid);
    qv_ref = interp1(z, qv_ref, zgrid);
    rho = interp1(z, rho, zgrid);
    w_9999 = interp1(masked.z, masked.w_9999, zgrid);
    p = interp1(z, p, zgrid) * 1e2;
    qcum = cumsum(dzqn_ref.*dz);
    qtilde = 1 - qcum/qcum(end);
    qtilde = interp1(z, qtilde, zgrid);
    w_99999 = interp1(masked.z, masked.w_9999, zgrid);
    z = zgrid;

    % Calculate saturation specific humidity profile
    qstar_ref = SOM_qsat(tabs_ref, p);

    % Find index of level of peak mass flux
    imax = find(mass_flux_cond9999 == max(mass_flux_cond9999), 1);
    i220 = find(tabs_ref <= 220, 1);

    % Create liquid/ice water static energy and total water arrays
    hL = NaN(1, numel(z));
    qt = NaN(1, numel(z));
    w2 = NaN(1, numel(z));
    w2_reda = NaN(1, numel(z));

    % Create environmental profiles
    hLe = SAM_c_p * tabs_ref + SAM_g * z;
    qte = RH * qstar_ref;

    % Find location of lower boundary (1 km)
    ilcl = find(z >= 1000, 1);
    
    % Iterate on temperature to find state that matches simulated
    % buoyancy profile at lower boundary
    Tzblow = tabs_ref(ilcl) - 10;
    Tzbhigh = tabs_ref(ilcl) + 10;
    qlow = SOM_qsat(Tzblow, p(ilcl));
    qhigh = SOM_qsat(Tzbhigh, p(ilcl));
    buoy_low = SAM_g * ((Tzblow - tabs_ref(ilcl)) ./ tabs_ref(ilcl) + ...
        0.608 * (qlow - qv_ref(ilcl)));
    buoy_high = SAM_g * ((Tzbhigh - tabs_ref(ilcl)) ./ tabs_ref(ilcl) + ...
        0.608 * (qhigh - qv_ref(ilcl)));
    buoy_target = masked_buoy9999(ilcl);
    buoy_anom(ii) = buoy_target;
    while (buoy_high - buoy_low > 0.0001)
        fprintf('%d\n', buoy_high - buoy_low);
        Tzbmid = 0.5 * (Tzblow + Tzbhigh);
        qmid = SOM_qsat(Tzbmid, p(ilcl));
        buoy_mid = SAM_g * ((Tzbmid - tabs_ref(ilcl)) ./ tabs_ref(ilcl) + ...
            0.608 * (qmid - qv_ref(ilcl)));
        if buoy_mid < buoy_target
            Tzblow = Tzbmid;
        else
            Tzbhigh = Tzbmid;
        end
        qlow = SOM_qsat(Tzblow, p(ilcl));
        qhigh = SOM_qsat(Tzbhigh, p(ilcl));
        buoy_low = SAM_g * ((Tzblow - tabs_ref(ilcl)) ./ tabs_ref(ilcl) + ...
            0.608 * (qlow - qv_ref(ilcl)));
        buoy_high = SAM_g * ((Tzbhigh - tabs_ref(ilcl)) ./ tabs_ref(ilcl) + ...
            0.608 * (qhigh - qv_ref(ilcl)));
    end

    Tzb = 0.5 * (Tzblow + Tzbhigh);
    qt(ilcl) = SOM_qsat(Tzb, p(ilcl));
    hL(ilcl) = SAM_c_p * Tzb + SAM_g * z(ilcl);
    w2(ilcl) = 0;
    w2_reda(ilcl) = 0;

    % Integrate
    for zz = (ilcl+1):numel(z)
        dhLdz = -eps * (hL(zz-1) - hLe(zz-1));
        hL(zz) = hL(zz-1) + (z(zz) - z(zz-1)) * dhLdz;
        dqtdz = -eps * (qt(zz-1) - qte(zz));
        qt(zz) = qt(zz-1) + (z(zz) - z(zz-1)) * dqtdz;
    end
    
    % Invert for thermodynamic state
    inv = SOM_getState(hL, qt, zeros(size(qt)), z, p);
    tabs = inv.T;
    qn = inv.qn;
    qv = qt - qn;

    % Calculate buoyancy, MSE, and saturation MSE
    buoy = SAM_g * ((tabs - tabs_ref) ./ tabs_ref + ...
        0.608 * (qv - qv_ref) - qn);

    % Integrate buoyancy
    for zz = (ilcl+1):numel(z)
        dw2dz = 2 * (buoy(zz-1) - 2*eps*w2(zz-1));
        w2(zz) = w2(zz-1) + (z(zz) - z(zz-1)) * dw2dz;
        dw2dz = 2 * (1/7 * buoy(zz-1)); %  - 2*eps*w2_reda(zz-1));
        w2_reda(zz) = w2_reda(zz-1) + (z(zz) - z(zz-1)) * dw2dz;
    end

    % Plot result
    figure(f1); hold on;
    l1 = [l1, plot(SAM_g * (tabs - tabs_ref) ./ tabs_ref, ...
        z/1e3, 'Color', cols(ii,:), 'LineWidth', 1.5)];
    plot(SAM_g * (tabs(imax) - tabs_ref(imax)) ./ tabs_ref(imax), ...
        z(imax)/1e3, 'o', 'MarkerFaceColor', cols(ii,:), ...
        'MarkerEdgeColor', cols(ii,:), 'MarkerSize', 10);
    plot(SAM_g * (tabs(i220) - tabs_ref(i220)) ./ tabs_ref(i220), ...
        z(i220)/1e3, 'x', 'MarkerFaceColor', cols(ii,:), ...
        'MarkerEdgeColor', cols(ii,:), 'MarkerSize', 10);
    
    figure(f2); hold on;
    l2 = [l2, plot(SAM_g * 0.608 * (qv - qv_ref), ...
        z/1e3, 'Color', cols(ii,:), 'LineWidth', 1.5)];
    plot(SAM_g * 0.608 * (qv(imax) - qv_ref(imax)), ...
        z(imax)/1e3, 'o', 'MarkerFaceColor', cols(ii,:), ...
        'MarkerEdgeColor', cols(ii,:), 'MarkerSize', 10);
    plot(SAM_g * 0.608 * (qv(i220) - qv_ref(i220)), ...
        z(i220)/1e3, 'x', 'MarkerFaceColor', cols(ii,:), ...
        'MarkerEdgeColor', cols(ii,:), 'MarkerSize', 10);

    figure(f3); hold on;
    l3 = [l3, plot(-SAM_g * qn, z/1e3, 'Color', cols(ii,:), 'LineWidth', 1.5)];
    plot(-SAM_g * qn(imax), z(imax)/1e3, 'o', ...
        'MarkerFaceColor', cols(ii,:), ...
        'MarkerEdgeColor', cols(ii,:), 'MarkerSize', 10);
    plot(-SAM_g * qn(i220), z(i220)/1e3, 'x', ...
        'MarkerFaceColor', cols(ii,:), ...
        'MarkerEdgeColor', cols(ii,:), 'MarkerSize', 10);
    
    figure(f4); hold on;
    l4 = [l4, plot(buoy, z/1e3, 'Color', cols(ii,:), 'LineWidth', 1.5)];
    plot(buoy(imax), z(imax)/1e3, 'o', ...
        'MarkerFaceColor', cols(ii,:), ...
        'MarkerEdgeColor', cols(ii,:), 'MarkerSize', 10);
    plot(buoy(i220), z(i220)/1e3, 'x', ...
        'MarkerFaceColor', cols(ii,:), ...
        'MarkerEdgeColor', cols(ii,:), 'MarkerSize', 10);
    
    figure(f5); hold on;
    l5 = [l5, plot(masked_buoy9999, z/1e3, ...
        'Color', cols(ii,:), 'LineWidth', 1.5)];
    plot(masked_buoy9999(imax), z(imax)/1e3, 'o', ...
        'MarkerFaceColor', cols(ii,:), ...
        'MarkerEdgeColor', cols(ii,:), 'MarkerSize', 10);
    plot(masked_buoy9999(i220), z(i220)/1e3, 'x', ...
        'MarkerFaceColor', cols(ii,:), ...
        'MarkerEdgeColor', cols(ii,:), 'MarkerSize', 10);
    
    % Plot pressure velocity composite members
    masked_omega = mass_flux_cond9999_all;
    masked_omega(masked_omega == 0) = nan;
    figure(standard_figure('half-page')); hold on;
    plot(SAM_g*masked_omega, p/1e2, 'Color', [0 0 0 0.02]);
    plot(SAM_g*mass_flux_cond9999, p/1e2, 'Color', 'black', 'LineWidth', 2);
    plot(rho(imax)*SAM_g*w_9999(imax), p(imax)/1e2, ...
        'ko', 'LineWidth', 2);
    xlabel('\omega (Pa s^{-1})');
    ylabel('Pressure (hPa)');
    ylim([0, 1000]);
    xlim([-10, 150]);
    set(gca, 'YDir', 'reverse');
    title(sprintf('Channel %d K', ss));
    pngandpdf(sprintf('FigureA4_%d', ss));

    % Compute buoyancy integrals
    buoyint_model_qtilde(ii) = trapz(z(ilcl:imax), buoy(ilcl:imax));
    buoyint_sim_qtilde(ii) = trapz(z(ilcl:imax), masked_buoy9999(ilcl:imax));
    buoyint_model_220K(ii) = trapz(z(ilcl:i220), buoy(ilcl:i220));
    buoyint_sim_220K(ii) = trapz(z(ilcl:i220), masked_buoy9999(ilcl:i220));
    
    % Compute vertical velocities
    Cd = 0.001;
    rho_qtilde(ii) = rho(imax);
    w_reda_qtilde(ii) = sqrt(w2_reda(imax));
    w_sticky_qtilde(ii) = sqrt(6/7 * buoy(imax) / Cd);
    w_cond_qtilde(ii) = mass_flux_cond9999(imax) / rho(imax);
    w_full_qtilde(ii) = sqrt(w2(imax));
    w_9999_qtilde(ii) = w_9999(imax);
    w_9999_max(ii) = max(w_9999);
    
    % Calculate vertical velocities in individual profiles
    % with peaks within 100 hPa of peak
    cmax = find(mass_flux_cond9999 == ...
       max(mass_flux_cond9999), 1);
    composite_pmax = p(cmax);
    w_cond_peak(ii) = 0;
    N = 0;
    tol = 75e2;
    for jj = 1:size(mass_flux_cond9999_all, 2)
       imax = find(mass_flux_cond9999_all(:,jj) == ...
           max(mass_flux_cond9999_all(:,jj)), 1);
       individual_pmax = p(imax);
       if abs(individual_pmax - composite_pmax) < tol
           w_cond_peak(ii) = w_cond_peak(ii) + ...
               mass_flux_cond9999_all(cmax,jj) / rho(cmax);
           N = N+1;
       end
    end
    w_cond_peak(ii) = w_cond_peak(ii) / N;
    w_cond_all(:,ii) = (mass_flux_cond9999_all(cmax,:) / rho(cmax))';

    % Save surface temperature perturbation
    Tzb_all(ii) = Tzb - tabs_ref(ilcl);
    
    % Iterate on temperature to find state that gives zero buoyancy
    % at cloud base
    Tzblow = tabs_ref(ilcl) - 10;
    Tzbhigh = tabs_ref(ilcl) + 10;
    qlow = SOM_qsat(Tzblow, p(ilcl));
    qhigh = SOM_qsat(Tzbhigh, p(ilcl));
    buoy_low = SAM_g * ((Tzblow - tabs_ref(ilcl)) ./ tabs_ref(ilcl) + ...
        0.608 * (qlow - qv_ref(ilcl)));
    buoy_high = SAM_g * ((Tzbhigh - tabs_ref(ilcl)) ./ tabs_ref(ilcl) + ...
        0.608 * (qhigh - qv_ref(ilcl)));
    buoy_target = 0;
    while (buoy_high - buoy_low > 0.0001)
        fprintf('%d\n', buoy_high - buoy_low);
        Tzbmid = 0.5 * (Tzblow + Tzbhigh);
        qmid = SOM_qsat(Tzbmid, p(ilcl));
        buoy_mid = SAM_g * ((Tzbmid - tabs_ref(ilcl)) ./ tabs_ref(ilcl) + ...
            0.608 * (qmid - qv_ref(ilcl)));
        if buoy_mid < buoy_target
            Tzblow = Tzbmid;
        else
            Tzbhigh = Tzbmid;
        end
        qlow = SOM_qsat(Tzblow, p(ilcl));
        qhigh = SOM_qsat(Tzbhigh, p(ilcl));
        buoy_low = SAM_g * ((Tzblow - tabs_ref(ilcl)) ./ tabs_ref(ilcl) + ...
            0.608 * (qlow - qv_ref(ilcl)));
        buoy_high = SAM_g * ((Tzbhigh - tabs_ref(ilcl)) ./ tabs_ref(ilcl) + ...
            0.608 * (qhigh - qv_ref(ilcl)));
    end
    
    % Compute buoyancy integrals with zero cloud-base buoyancy
    Tzb = 0.5 * (Tzblow + Tzbhigh);
    qt(ilcl) = SOM_qsat(Tzb, p(ilcl));
    hL(ilcl) = SAM_c_p * Tzb + SAM_g * z(ilcl);
    w2(ilcl) = 0;
    w2_reda(ilcl) = 0;

    % Integrate
    for zz = (ilcl+1):numel(z)
        dhLdz = -eps * (hL(zz-1) - hLe(zz-1));
        hL(zz) = hL(zz-1) + (z(zz) - z(zz-1)) * dhLdz;
        dqtdz = -eps * (qt(zz-1) - qte(zz));
        qt(zz) = qt(zz-1) + (z(zz) - z(zz-1)) * dqtdz;
    end
    
    % Invert for thermodynamic state
    inv = SOM_getState(hL, qt, zeros(size(qt)), z, p);
    tabs = inv.T;
    qn = inv.qn;
    qv = qt - qn;

    % Calculate buoyancy, MSE, and saturation MSE
    buoy = SAM_g * ((tabs - tabs_ref) ./ tabs_ref + ...
        0.608 * (qv - qv_ref) - qn);
    
    % Plot profile
    figure(f6); hold on;
    plot(buoy, z/1e3, 'Color', cols(ii,:), 'LineWidth', 1.5);
    plot(buoy(imax), z(imax)/1e3, 'o', ...
        'MarkerFaceColor', cols(ii,:), ...
        'MarkerEdgeColor', cols(ii,:), 'MarkerSize', 10);
    plot(buoy(i220), z(i220)/1e3, 'x', ...
        'MarkerFaceColor', cols(ii,:), ...
        'MarkerEdgeColor', cols(ii,:), 'MarkerSize', 10);
    
    lbl{ii} = sprintf('%d K', ss);
    ii = ii+1;
end

%%
figure(f1);
xlabel('Temperature perturbation (m s^{-2})');
ylabel('Height (km)');
ylim([0 15]);
xlim([0 0.2]);
legend(l1, lbl);

figure(f2);
xlabel('Virtual effect (m s^{-2})');
ylabel('Height (km)');
ylim([0 15]);
xlim([0 0.1]);
legend(l2, lbl);

figure(f3);
xlabel('Condensate loading (m s^{-2})');
ylabel('Height (km)');
ylim([0 15]);
xlim([-0.1 0]);
legend(l3, lbl);

figure(f4);
xlabel('Plume model buoyancy (m s^{-2})');
ylabel('Height (km)');
ylim([0 15]);
xlim([-0.02 0.15]);
title('(b)');
text(0.06, 2.5, sprintf('Diagnosed cloud-base\nbuoyancy'));
plot([0 0], ylim, 'k--');
pngandpdf('Figure9B');

figure(f5);
xlabel('Simulated buoyancy (m s^{-2})');
ylabel('Height (km)');
ylim([0 15]);
xlim([0 0.15]);
leg = legend(l5, lbl, 'Location', 'southeast');
leg.Box = 'off';
leg.ItemTokenSize = [10,1];
title('(a)');
pngandpdf('Figure9A');

figure(f6);
xlabel('Plume model buoyancy (m s^{-2})');
ylabel('Height (km)');
ylim([0 15]);
xlim([-0.02 0.15]);
title('(d)');
text(0.07, 2.5, sprintf('Zero cloud-base\nbuoyancy'));
plot([0 0], ylim, 'k--');
pngandpdf('Figure9D');

%%
figure(standard_figure('half-page')); hold on;
plot(sst, buoyint_model_qtilde, 'k-o', 'LineWidth', 1.5);
plot(sst, buoyint_sim_qtilde, 'r-o', 'LineWidth', 1.5);
xlabel('SST (K)');
ylabel('Buoyancy integral (J kg^{-1})');
leg = legend('Plume model', 'Simulation 99.99th p''tile', ...
    'Location', 'northwest');
leg.Box = 'off';
leg.ItemTokenSize = [10,1];
title('(a)');
yl = ylim;
ylim([0 yl(2)]);
pngandpdf('Figure8A');

figure(standard_figure('half-page')); hold on;
plot(sst, buoyint_model_220K, 'k-o', 'LineWidth', 1.5);
plot(sst, buoyint_sim_220K, 'r-o', 'LineWidth', 1.5);
xlabel('SST (K)');
ylabel('Buoyancy integral (J kg^{-1})');
legend('Plume model', 'Simulation 99.99th percentile', ...
    'Location', 'northwest');

figure(standard_figure('half-page')); hold on;
plot(sst, w_full_qtilde, 'k-o', 'LineWidth', 1.5);
plot(sst, w_cond_qtilde, '-o', 'LineWidth', 1.5, 'Color', [0 0.5 0]);
plot(sst, w_9999_max, 'b-o', 'LineWidth', 1.5);
xlabel('SST (K)');
ylabel('Vertical velocity (m s^{-1})');
legend('Plume model, w^2 equation', ...
    'Simulation conditional average', 'Simulation 99.99th percentile (max)', ...
    'Location', 'northwest');

figure(standard_figure('half-page')); hold on;
plot(sst, w_full_qtilde, 'k-o', 'LineWidth', 1.5);
plot(sst, w_cond_qtilde, '-o', 'LineWidth', 1.5, 'Color', [0 0.5 0]);
plot(sst, w_9999_qtilde, 'b-o', 'LineWidth', 1.5);
xlabel('SST (K)');
ylabel('Vertical velocity (m s^{-1})');
legend('Plume model, w^2 equation', ...
    'Simulation conditional average', 'Simulation 99.99th percentile w', ...
    'Location', 'northwest');

figure(standard_figure('half-page')); hold on;
plot(sst, w_reda_qtilde, 'k-o', 'LineWidth', 1.5);
plot(sst, w_cond_qtilde, '-o', 'LineWidth', 1.5, 'Color', [0 0.5 0]);
plot(sst, w_9999_qtilde, 'b-o', 'LineWidth', 1.5);
xlabel('SST (K)');
ylabel('Vertical velocity (m s^{-1})');
leg = legend('Plume model, w^2/(2BI) = 0.14', ...
    'Simulation conditional average', ...
    'Simulation 99.99th p''tile', ...
    'Location', 'northwest');
leg.Box = 'off';
leg.ItemTokenSize = [10,1];
title('(b)');
yl = ylim;
ylim([0 yl(2)]);
pngandpdf('Figure8');

figure(standard_figure('half-page')); hold on;
plot(sst, w_reda_qtilde, 'k-o', 'LineWidth', 1.5);
plot(sst, w_cond_qtilde, '-o', 'LineWidth', 1.5, 'Color', [0 0.5 0]);
plot(sst, w_cond_peak, '--o', 'LineWidth', 1.5, 'Color', [0 0.5 0]);
plot(sst, w_9999_qtilde, 'b-o', 'LineWidth', 1.5);
xlabel('SST (K)');
ylabel('Vertical velocity (m s^{-1})');
leg = legend('Plume model, w^2/(2BI) = 0.14', ...
    'Simulation conditional average', ...
    'Subsampled conditional average', ...
    'Simulation 99.99th p''tile', ...
    'Location', 'northwest');
leg.Box = 'off';
leg.ItemTokenSize = [10,1];
title('(b)');
yl = ylim;
ylim([0 yl(2)]);
pngandpdf('Figure8B_app');

figure(standard_figure('half-page')); hold on;
plot(sst, -rho_qtilde.*w_reda_qtilde/9.81, 'k-o', 'LineWidth', 1.5);
plot(sst, -rho_qtilde.*w_cond_qtilde/9.81, '-o', 'LineWidth', 1.5, 'Color', [0 0.5 0]);
plot(sst, -rho_qtilde.*w_9999_qtilde/9.81, 'b-o', 'LineWidth', 1.5);
xlabel('SST (K)');
ylabel('Pressure velocity (Pa s^{-1})');
leg = legend('Plume model, \eta = 1/7', ...
    'Simulation conditional average', ...
    'Simulation 99.99th p''tile', ...
    'Location', 'northwest');
leg.Box = 'off';
leg.ItemTokenSize = [10,1];
set(gca, 'YDir', 'reverse');

figure(standard_figure('half-page')); hold on;
plot(sst, w_reda_qtilde, 'k-o', 'LineWidth', 1.5);
plot(sst, w_cond_qtilde, '-o', 'LineWidth', 1.5, 'Color', [0 0.5 0]);
xlabel('SST (K)');
ylabel('Vertical velocity (m s^{-1})');
legend('Plume model, w^2 equation with \eta = 1/7', ...
    'Simulation conditional average', ...
    'Location', 'northwest');

figure(standard_figure('half-page')); hold on;
plot(sst, w_sticky_qtilde, 'k-o', 'LineWidth', 1.5);
plot(sst, w_cond_qtilde, '-o', 'LineWidth', 1.5, 'Color', [0 0.5 0]);
plot(sst, w_9999_qtilde, 'b--o', 'LineWidth', 1.5);
xlabel('SST (K)');
ylabel('Vertical velocity (m s^{-1})');
legend('Plume model, sticky limit', ...
    'Simulation conditional average', 'Simulation 99.99th percentile w at p_{max}', ...
    'Location', 'northwest');

figure(standard_figure('half-page')); hold on;
plot(sst, buoy_anom, 'k-o', 'LineWidth', 2);
xlabel('SST (K)');
ylabel('Cloud base buoyancy (m s^{-2})');
title('(c)');
pngandpdf('Figure9C');

figure(standard_figure('half-page')); hold on;
ii = 1;
for ss = sst
    [N, edges] = histcounts(w_cond_all(:,ii)/w_9999_qtilde(ii), ...
        'Normalization', 'pdf');
    centers = 0.5 * (edges(2:end) + edges(1:end-1));
    plot(centers, N, 'LineWidth', 1.5, 'Color', cols(ii,:));
    ii = ii+1;
end
xlabel('w / w_{99.99} (nondim)');
ylabel('Frequency (nondim)');
leg = legend(gca, lbl, 'Units', 'normalized', 'Location', 'northwest');
leg.ItemTokenSize = [10,1];
leg.Box = 'off';
title('(a) Channel');
pngandpdf('Figure9_A1');
