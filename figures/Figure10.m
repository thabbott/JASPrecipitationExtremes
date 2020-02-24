%% ABC2020 Figure 10

% Setup
addpath(genpath('..\matlab'));
latex = 0;
set_plot_styling;
SAM_defineConstants;
SOM_defineConstants;

%% Vary several parameters
% Add points for dynamic mode with cloud base MSE anomalies 
% diagnosed from simulations
% (and maybe error bars derived based on uncertainty in other parameters?)

% entrainment rates
params = [0.5e-3, 0, 0.9;
          0.5e-3, 0.01, 0.9;
          0.5e-3, 0.05, 0.9;
          0.5e-3, 0.1, 0.9];
      
lstyle = {'-', '-', '-', '-'};
cols = [0 0 0;
    0 0 1;
    0 0.5 0;
    1 0 0];
lbl = {'B_b = 0', 'B_b = 0.01 m s^{-2}', 'B_b = 0.05 m s^{-2}', ...
    'B_b = 0.1 m s^{-2}'};
          

sRH = 0.8;
Lv = SAM_L_c; % J kg-1
Rv = SAM_R_v; % J kg-1 K-1
Ra = SAM_R; % J kg-1 K-1
cp = SAM_c_p; % J kg-1 K-1
g = SAM_g; % N kg-1
fm3 = standard_figure('half-page'); hold on;
fm2 = standard_figure('half-page'); hold on;
fm1 = standard_figure('half-page'); hold on;
f0 = standard_figure('half-page'); hold on;
f1 = standard_figure('half-page'); hold on;
f2 = standard_figure('half-page'); hold on;
l = [];
l1 = [];
l2 = [];

for ii = 1:size(params, 1)
    
    % Initialize surface properties and pressure grid
    Ts = linspace(270, 330, 13);
    ps = 1e5*ones(size(Ts));
    qvs = sRH*SOM_qsat(Ts, ps);
    p = linspace(1e5,100,1000);
    p = reshape(p, [1 1 numel(p)]) .* ones(size(ps));

    % Loop over surface temperatures, calculate buoyancy profiles and integrals
    epss = params(ii,1);
    B_b = params(ii,2);
    RH = params(ii,3);
    epsw = 0.15e-3;
    deps = epss - epsw;
    qtilde_max = 0.35; % nondim
    Lv = SAM_L_c; % J kg-1
    Rv = SAM_R_v; % J kg-1 K-1
    cp = SAM_c_p; % J kg-1 K-1
    g = SAM_g; % N kg-1
    Hq = 3e3; % m
    dT = zeros(size(p));
    BI = zeros(size(p,1), size(p,2));
    rho_qt = zeros(size(p,1), size(p,2));
    T_qt = zeros(size(p,1), size(p,2));
    p_qt = zeros(size(p,1), size(p,2));

	% Calculate environmental profiles
    fprintf('Doing zero-buoyancy plume model calculations\n');
    [T, ~, qstar, p, z] = CC_zbplume(Ts, qvs, p, 1, ...
        RH, epss);
    fprintf('Done\n');
    rho = p./(SAM_R * T);
    
    % Calculate surface temperature
    Tsurf = Ts;
    
    % Calculate surface saturation specific humidity
    qstar_surf = SOM_qsat(Tsurf, ps);

    % Calculate qtilde
    qtilde = qstar ./ qstar_surf;
    
    for j = 1:size(p,2)
        for i = 1:size(p,1)
            
            % Find qtilde = 0.35
            iqt = find(qtilde(i,j,:) <= qtilde_max, 1, 'first');

            for k = 2:iqt

                % Calculate exponential kernel
                kernel = exp(-epsw * (z(i,j,k) - squeeze(z(i,j,1:k))));
                cloud = (z(i,j,1:k) >= 1e3);
                kernel(~cloud) = 0;
                if any(cloud)
                    ilcl = find(cloud, 1);
                    weight = exp(-epsw * (z(i,j,k) - z(i,j,ilcl)));
                    Tb = T(i,j,ilcl);
                    dTb = B_b * Tb / SAM_g;
                    dhb = cp * dTb + ...
                        Lv * (SOM_qsat(Tb + dTb, p(i,j,ilcl)) - ...
                        SOM_qsat(Tb, p(i,j,ilcl)));
                    disp(z(i,j,ilcl));
                else
                    dhb = 0;
                    weight = 0;
                end

                % Calculate temperature perturbation
                dT(i,j,k) = 1 ./ ...
                    (1 + Lv^2 * qstar(i,j,k) ./ (cp * Rv * T(i,j,k).^2)) .* ...
                    (weight * dhb/cp + ...
                    deps * (1 - RH) * Lv / cp * ...
                    trapz(squeeze(z(i,j,1:k)), kernel.*squeeze(qstar(i,j,1:k))));

            end

            % Calculate buoyancy integral
            BI(i,j) = trapz(squeeze(z(i,j,1:iqt)), ...
                g*squeeze(dT(i,j,1:iqt)./T(i,j,1:iqt)));
            rho_qt(i,j) = rho(i,j,iqt);
            T_qt(i,j) = T(i,j,iqt);
            p_qt(i,j) = p(i,j,iqt);
        end
    end
    
    figure(fm3); hold on;
    plot(Tsurf, p_qt, lstyle{ii}, 'LineWidth', 1.5, ...
        'Color', cols(ii,:));
    
    figure(fm2); hold on;
    plot(Tsurf, rho_qt, lstyle{ii}, 'LineWidth', 1.5, ...
        'Color', cols(ii,:));

    figure(fm1); hold on;
    l1 = [l1, ...
        plot(Tsurf, sqrt(2*0.14*BI), lstyle{ii}, 'LineWidth', 1.5, ...
        'Color', cols(ii,:))];
    
    omega = -g*rho_qt.*sqrt(2*0.14*BI);
    figure(f0); hold on;
    l2 = [l2, ...
        plot(Tsurf, omega, lstyle{ii}, 'LineWidth', 1.5, ...
        'Color', cols(ii,:))];
    Tsm = 0.5*(Tsurf(1:end-1) + Tsurf(2:end));
    scal = @(x) 200 * (x(2:end) - x(1:end-1)) ./ (x(2:end) + x(1:end-1)) ./ ...
    diff(Tsurf);
    p = polyfit(Tsurf, log(omega), 4);
    pd = polyder(p);   
    plot(Tsurf, exp(polyval(p, Tsurf)), [lstyle{ii}, '-'], ...
        'LineWidth', 0.5, 'Color', cols(ii,:));
    
    figure(f1); hold on;
    l = [l, ...
        plot(Tsurf, 100*polyval(pd, Tsurf), lstyle{ii}, ...
        'LineWidth', 1.5, 'Color', cols(ii,:))];
    
    figure(f2); hold on;
    plot(Tsm, scal(sqrt(2.0*1.4*BI)), lstyle{ii}, ...
        'LineWidth', 1.5, 'Color', cols(ii,:));

end

figure(fm3); hold on;
xlabel('Surface air T (K)');
ylabel('Pressure (Pa)');
xlim([min(Tsurf), max(Tsurf)]);
xlim([275, 325]);

figure(fm2); hold on;
xlabel('Surface air T (K)');
ylabel('Density (kg m^{-3})');
xlim([min(Tsurf), max(Tsurf)]);
xlim([275, 325]);

figure(fm1); hold on;
xlabel('Surface air T (K)');
ylabel('Vertical velocity w (m s^{-1})');
xlim([min(Tsurf), max(Tsurf)]);
xlim([275, 325]);

figure(f0); hold on;
xlabel('Surface air T (K)');
ylabel('Pressure velocity \omega (Pa s^{-1})');
xlim([min(Tsurf), max(Tsurf)]);
set(gca, 'YDir', 'reverse');
xlim([275, 325]);

figure(f1); hold on;
plot(Tsm, zeros(size(Tsm)), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
xlabel('Surface air T (K)');
ylabel('Dynamic mode (% K^{-1})');
xlim([min(Tsurf), max(Tsurf)]);
xlim([275, 325]);
title('(d)');

figure(f2); hold on;
plot(Tsm, zeros(size(Tsm)), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
xlabel('Surface air T (K)');
ylabel('Change in w (% K^{-1})');
xlim([275, 325]);

%% Add dots from simulations
sst = 280:5:305;
get_export_fname = @(sst, fname) sprintf('../data/ch_cam%dri0/%s', ...
    sst, fname);
ii = 1;
nboot = 1000;
Tsurf = zeros(numel(sst), 1);
peak_mass_flux = zeros(numel(sst), 1);
peak_w = zeros(numel(sst), 1);
condm = zeros(numel(sst)-1,1);
condm_ref = zeros(numel(sst)-1,1);
condm_boot = zeros(numel(sst)-1,nboot);
condm_boot_ref = zeros(numel(sst)-1,nboot);
for ss = sst
    load(get_export_fname(ss, 'extremes9999.mat'), 'mass_flux_cond9999', ...
        'precip9999', 'p', 'rho', 'dz', 'dzqn_ref', 'tabs_ref');
    Tsurf(ii) = tabs_ref(1);
    peak_mass_flux(ii) = max(mass_flux_cond9999);
    ipeak = find(mass_flux_cond9999 == peak_mass_flux(ii), 1);
    peak_w(ii) = peak_mass_flux(ii) / rho(ipeak);
    load(get_export_fname(ss, 'extremes9999_boot.mat'), ...
        'mass_flux_cond9999_boot');
    if ii < numel(sst)        
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

figure(fm1); hold on;
plot(Tsurf, peak_w, '*', 'LineWidth', 1.5, 'Color', [0 0.5 0]);
leg = legend(l1, lbl, 'Position', [0.13, 0.5, 0.3380, 0.4391]);
leg.ItemTokenSize = [10,1];
leg.Box = 'off';
plot(300, 2, '*', 'Color', [0 0.5 0], 'LineWidth', 1.5);
text(302, 2, 'Simulations', 'Color', [0 0.5 0]);

figure(f0); hold on;
l2 = [l2, ...
    plot(Tsurf, -SAM_g * peak_mass_flux, '*', 'LineWidth', 1.5, ...
    'Color', [0 0.5 0])];

figure(f1); hold on;
scal = @(x,y) 2*100*(x - y)./(x + y)./diff(Tsurf);
Tsm = 0.5*(Tsurf(1:end-1) + Tsurf(2:end));
condm_boot = scal(condm_boot, condm_boot_ref(1:end,:));
condm_boot_ext = [min(condm_boot, [], 2) max(condm_boot, [], 2)];
condm_boot = std(condm_boot, 0, 2);
errorbar(Tsm, scal(condm(1:end), condm_ref(1:end)), ...
    scal(condm(1:end), condm_ref(1:end)) - condm_boot_ext(:,1),...
    -(scal(condm(1:end), condm_ref(1:end)) - condm_boot_ext(:,2)), ...
    '*', ...
    'Color', [0 0.5 0], 'LineWidth', 1.5);
plot(280, -1.1, '*', 'Color', [0 0.5 0], 'LineWidth', 1.5);
text(282, -1.1, 'Simulations', 'Color', [0 0.5 0]);
leg = legend(l, lbl, 'Location', 'northeast');
leg.ItemTokenSize = [20,1];
leg.Box = 'off';

figure(f2); hold on;
plot(Tsm, zeros(size(Tsm)), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
xlabel('Surface air T (K)');
ylabel('Change in w (% K^{-1})');
xlim([min(Tsurf), max(Tsurf)]);
leg = legend(l, lbl, 'Location', 'northeast');
leg.ItemTokenSize = [10,1];
leg.Box = 'off';
xlim([275, 325]);
ylim([-2, 8]);

%% Save figures
figure(fm1); hold on;
title('(a)');
pngandpdf('Figure10A');

figure(f0); hold on;
title('(b)');
pngandpdf('Figure10B');

figure(f1); hold on;
title('(c)');
pngandpdf('Figure10C');

%% Re-do calculations with varying weak entrainment rates

% entrainment rates
params = [0.5e-3, 0.05, 0.9, 0;
          0.5e-3, 0.05, 0.9, 0.1e-3;
          0.5e-3, 0.05, 0.9, 0.2e-3;
          0.5e-3, 0.05, 0.9, 0.3e-3];
      
lstyle = {'-', '-', '-', '-'};
cols = copper(size(params, 1));
lbl = {'\epsilon_w = 0 km^{-1}', '\epsilon_w = 0.1 km^{-1}', ...
    '\epsilon_w = 0.2 km^{-1}', '\epsilon_w = 0.3 km^{-1}'};
          

sRH = 0.8;
Lv = SAM_L_c; % J kg-1
Rv = SAM_R_v; % J kg-1 K-1
Ra = SAM_R; % J kg-1 K-1
cp = SAM_c_p; % J kg-1 K-1
g = SAM_g; % N kg-1
fm3 = standard_figure('half-page'); hold on;
fm2 = standard_figure('half-page'); hold on;
fm1 = standard_figure('half-page'); hold on;
f0 = standard_figure('half-page'); hold on;
f1 = standard_figure('half-page'); hold on;
f2 = standard_figure('half-page'); hold on;
l = [];
l2 = [];

for ii = 1:size(params, 1)
    
    % Initialize surface properties and pressure grid
    Ts = linspace(270, 330, 13);
    ps = 1e5*ones(size(Ts));
    qvs = sRH*SOM_qsat(Ts, ps);
    p = linspace(1e5,100,1000);
    p = reshape(p, [1 1 numel(p)]) .* ones(size(ps));

    % Loop over surface temperatures, calculate buoyancy profiles and integrals
    epss = params(ii,1);
    B_b = params(ii,2);
    RH = params(ii,3);
    epsw = params(ii,4);
    deps = epss - epsw;
    qtilde_max = 0.35; % nondim
    Lv = SAM_L_c; % J kg-1
    Rv = SAM_R_v; % J kg-1 K-1
    cp = SAM_c_p; % J kg-1 K-1
    g = SAM_g; % N kg-1
    Hq = 3e3; % m
    dT = zeros(size(p));
    BI = zeros(size(p,1), size(p,2));
    rho_qt = zeros(size(p,1), size(p,2));
    T_qt = zeros(size(p,1), size(p,2));
    p_qt = zeros(size(p,1), size(p,2));

	% Calculate environmental profiles
    fprintf('Doing zero-buoyancy plume model calculations\n');
    [T, ~, qstar, p, z] = CC_zbplume(Ts, qvs, p, 1, ...
        RH, epss);
    fprintf('Done\n');
    rho = p./(SAM_R * T);
    
    % Calculate surface temperature
    Tsurf = Ts;
    
    % Calculate surface saturation specific humidity
    qstar_surf = SOM_qsat(Tsurf, ps);

    % Calculate qtilde
    qtilde = qstar ./ qstar_surf;
    
    for j = 1:size(p,2)
        for i = 1:size(p,1)
            
            % Find qtilde = 0.35
            iqt = find(qtilde(i,j,:) <= qtilde_max, 1, 'first');

            for k = 2:iqt

                % Calculate exponential kernel
                kernel = exp(-epsw * (z(i,j,k) - squeeze(z(i,j,1:k))));
                cloud = (z(i,j,1:k) >= 1e3);
                kernel(~cloud) = 0;
                if any(cloud)
                    ilcl = find(cloud, 1);
                    weight = exp(-epsw * (z(i,j,k) - z(i,j,ilcl)));
                    Tb = T(i,j,ilcl);
                    dTb = B_b * Tb / SAM_g;
                    dhb = cp * dTb + ...
                        Lv * (SOM_qsat(Tb + dTb, p(i,j,ilcl)) - ...
                        SOM_qsat(Tb, p(i,j,ilcl)));
                    disp(z(i,j,ilcl));
                else
                    dhb = 0;
                    weight = 0;
                end

                % Calculate temperature perturbation
                dT(i,j,k) = 1 ./ ...
                    (1 + Lv^2 * qstar(i,j,k) ./ (cp * Rv * T(i,j,k).^2)) .* ...
                    (weight * dhb/cp + ...
                    deps * (1 - RH) * Lv / cp * ...
                    trapz(squeeze(z(i,j,1:k)), kernel.*squeeze(qstar(i,j,1:k))));

            end

            % Calculate buoyancy integral
            BI(i,j) = trapz(squeeze(z(i,j,1:iqt)), ...
                g*squeeze(dT(i,j,1:iqt)./T(i,j,1:iqt)));
            rho_qt(i,j) = rho(i,j,iqt);
            T_qt(i,j) = T(i,j,iqt);
            p_qt(i,j) = p(i,j,iqt);
        end
    end
    
    figure(fm3); hold on;
    plot(Tsurf, p_qt, lstyle{ii}, 'LineWidth', 1.5, ...
        'Color', cols(ii,:));
    
    figure(fm2); hold on;
    plot(Tsurf, rho_qt, lstyle{ii}, 'LineWidth', 1.5, ...
        'Color', cols(ii,:));

    figure(fm1); hold on;
    plot(Tsurf, sqrt(2*0.14*BI), lstyle{ii}, 'LineWidth', 1.5, ...
        'Color', cols(ii,:));
    
    omega = -g*rho_qt.*sqrt(2*0.14*BI);
    figure(f0); hold on;
    l2 = [l2, ...
        plot(Tsurf, omega, lstyle{ii}, 'LineWidth', 1.5, ...
        'Color', cols(ii,:))];
    Tsm = 0.5*(Tsurf(1:end-1) + Tsurf(2:end));
    scal = @(x) 200 * (x(2:end) - x(1:end-1)) ./ (x(2:end) + x(1:end-1)) ./ ...
    diff(Tsurf);
    p = polyfit(Tsurf, log(omega), 4);
    pd = polyder(p);   
    plot(Tsurf, exp(polyval(p, Tsurf)), [lstyle{ii}, '-'], ...
        'LineWidth', 0.5, 'Color', cols(ii,:));
    
    figure(f1); hold on;
    l = [l, ...
        plot(Tsurf, 100*polyval(pd, Tsurf), lstyle{ii}, ...
        'LineWidth', 1.5, 'Color', cols(ii,:))];
    
    figure(f2); hold on;
    plot(Tsm, scal(sqrt(2.0*1.4*BI)), lstyle{ii}, ...
        'LineWidth', 1.5, 'Color', cols(ii,:));

end

figure(fm3); hold on;
xlabel('Surface air T (K)');
ylabel('Pressure (Pa)');
xlim([min(Tsurf), max(Tsurf)]);
xlim([275, 325]);

figure(fm2); hold on;
xlabel('Surface air T (K)');
ylabel('Density (kg m^{-3})');
xlim([min(Tsurf), max(Tsurf)]);
xlim([275, 325]);

figure(fm1); hold on;
xlabel('Surface air T (K)');
ylabel('Vertical velocity (m s^{-1})');
xlim([min(Tsurf), max(Tsurf)]);
xlim([275, 325]);

figure(f0); hold on;
xlabel('Surface air T (K)');
ylabel('Pressure velocity (Pa s^{-1})');
xlim([min(Tsurf), max(Tsurf)]);
set(gca, 'YDir', 'reverse');
xlim([275, 325]);

figure(f1); hold on;
plot(Tsm, zeros(size(Tsm)), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
xlabel('Surface air T (K)');
ylabel('Dynamic mode (% K^{-1})');
xlim([min(Tsurf), max(Tsurf)]);
xlim([275, 325]);
title('(d)');

figure(f2); hold on;
plot(Tsm, zeros(size(Tsm)), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
xlabel('Surface air T (K)');
ylabel('Change in w (% K^{-1})');
xlim([275, 325]);

%% Add dots from simulations
sst = 280:5:305;
tas = [];
get_export_fname = @(sst, fname) sprintf('../data/ch_cam%dri0/%s', ...
    sst, fname);
ii = 1;
nboot = 1000;
Tsurf = zeros(numel(sst), 1);
peak_mass_flux = zeros(numel(sst), 1);
peak_w = zeros(numel(sst), 1);
condm = zeros(numel(sst)-1,1);
condm_ref = zeros(numel(sst)-1,1);
condm_boot = zeros(numel(sst)-1,nboot);
condm_boot_ref = zeros(numel(sst)-1,nboot);
for ss = sst
    load(get_export_fname(ss, 'extremes9999.mat'), 'mass_flux_cond9999', ...
        'precip9999', 'p', 'rho', 'dz', 'dzqn_ref', 'tabs_ref');
    Tsurf(ii) = tabs_ref(1);
    peak_mass_flux(ii) = max(mass_flux_cond9999);
    ipeak = find(mass_flux_cond9999 == peak_mass_flux(ii), 1);
    peak_w(ii) = peak_mass_flux(ii) / rho(ipeak);
    load(get_export_fname(ss, 'extremes9999_boot.mat'), ...
        'mass_flux_cond9999_boot');
    if ii < numel(sst)        
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

figure(fm1); hold on;
plot(Tsurf, peak_w, '*', 'LineWidth', 1.5, 'Color', [0 0.5 0]);

figure(f0); hold on;
l2 = [l2, ...
    plot(Tsurf, -SAM_g * peak_mass_flux, '*', 'LineWidth', 1.5, ...
    'Color', [0 0.5 0])];

figure(f1); hold on;
scal = @(x,y) 2*100*(x - y)./(x + y)./diff(Tsurf);
Tsm = 0.5*(Tsurf(1:end-1) + Tsurf(2:end));
condm_boot = scal(condm_boot, condm_boot_ref(1:end,:));
condm_boot_ext = [min(condm_boot, [], 2) max(condm_boot, [], 2)];
condm_boot = std(condm_boot, 0, 2);
errorbar(Tsm, scal(condm(1:end), condm_ref(1:end)), ...
    scal(condm(1:end), condm_ref(1:end)) - condm_boot_ext(:,1),...
    -(scal(condm(1:end), condm_ref(1:end)) - condm_boot_ext(:,2)), ...
    '*', ...
    'Color', [0 0.5 0], 'LineWidth', 1.5);
plot(280, -1.1, '*', 'Color', [0 0.5 0], 'LineWidth', 1.5);
text(282, -1.1, 'Simulations', 'Color', [0 0.5 0]);
leg = legend(l, lbl, 'Location', 'northeast');
leg.ItemTokenSize = [20,1];
leg.Box = 'off';
ylim([-2, 8]);

figure(f2); hold on;
plot(Tsm, zeros(size(Tsm)), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
xlabel('Surface air T (K)');
ylabel('Change in w (% K^{-1})');
xlim([min(Tsurf), max(Tsurf)]);
leg = legend(l, lbl, 'Location', 'northeast');
leg.ItemTokenSize = [10,1];
leg.Box = 'off';
xlim([275, 325]);

%% Save figures
figure(f1);
title('(d)');
pngandpdf('Figure10D');
