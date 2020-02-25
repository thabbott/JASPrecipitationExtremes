%% ABC2020 Figure A5

%% Setup
addpath(genpath('../matlab'));
latex = 0;
set_plot_styling;
SAM_defineConstants;
SOM_defineConstants;

% Calculate metrics for channel simulations
sst = 280:5:305;        % SST
fnamet = ['../Data/ch_cam%dri0/%s'];
zgrid = 50:50:20e3;          % Vertical grid (m)
ii = 1;
err_channel = zeros(numel(sst),1);
for ss = sst
    
    % Load profiles
    load(sprintf(fnamet, ss, 'extremes9999.mat'), ...
        'mass_flux_cond9999', 'z', 'p', ...
        'mass_flux_cond9999_all');
    
    % Interpolate onto high-resolution grid
    mass_flux_cond9999 = interp1(z, mass_flux_cond9999, zgrid)';
    mass_flux_cond9999_all = interp1(z, mass_flux_cond9999_all, zgrid);
    p = interp1(z, p, zgrid) * 1e2;
    
    % Calculate metric
    err = -trapz(p, ...
            abs(mass_flux_cond9999 - mass_flux_cond9999_all)) ./ sqrt(...
            trapz(p, mass_flux_cond9999) .* ...
            trapz(p, mass_flux_cond9999_all));
    err_channel(ii) = mean(err);
    
    ii = ii+1;
end

% Calculate metrics for small-domain simulations
sst = 280:5:305;        % SST
fnamet = ['../Data/RCEMIP_SST%d/%s'];
zgrid = 50:50:20e3;          % Vertical grid (m)
ii = 1;
err_small = zeros(numel(sst),1);
for ss = sst
    
    % Load profiles
    load(sprintf(fnamet, ss, 'extremes9999.mat'), ...
        'mass_flux_cond9999', 'z', 'p', ...
        'mass_flux_cond9999_all');
    
    % Interpolate onto high-resolution grid
    mass_flux_cond9999 = interp1(z, mass_flux_cond9999, zgrid)';
    mass_flux_cond9999_all = interp1(z, mass_flux_cond9999_all, zgrid);
    p = interp1(z, p, zgrid) * 1e2;
    
    % Calculate metric
    err = -trapz(p, ...
            abs(mass_flux_cond9999 - mass_flux_cond9999_all)) ./ sqrt(...
            trapz(p, mass_flux_cond9999) .* ...
            trapz(p, mass_flux_cond9999_all));
    err_small(ii) = mean(err);
    
    ii = ii+1;
end

% Plot
figure(standard_figure('half-page')); hold on;
plot(sst, err_channel, 'k', 'LineWidth', 1.5);
plot(sst, err_small, 'k--', 'LineWidth', 1.5);
xlabel('SST (K)');
ylabel('Mean collapse error (nondim.)');
leg = legend('Channel', 'Small-domain', 'Location', 'southeast');
leg.Box = 'off';
ylim([0, 1]);
pngandpdf('FigureA5');
