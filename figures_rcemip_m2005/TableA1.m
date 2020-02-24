%% ABC2020 Table 1

% Setup
addpath(genpath('..\matlab'));
SAM_defineConstants;
SOM_defineConstants;

sst = [280 285 290 295 300 305];
cols = viridis(numel(sst));
sstm = 0.5*(sst(2:end) + sst(1:end-1));
get_export_fname = @(sst, fname) sprintf('../data/RCEMIP_SST%d_M2005/%s', ...
    sst, fname);

edges_p = Inf * ones(numel(sst), numel(sst));
edges_q = Inf * ones(numel(sst), numel(sst));
jj = 1;
qgrid = linspace(1,0,1000);
pgrid = linspace(1e5,0,1000);
for ss = sst
    ii = 1;
    for ss2 = sst
        if ss2 == ss 
            ii = ii+1;
            continue; 
        end
        
        dat = load(get_export_fname(ss, 'extremes9999.mat'));
        dat2 = load(get_export_fname(ss2, 'extremes9999.mat'));
        dat.p = 1e2*dat.p;
        dat2.p = 1e2*dat2.p;
        
        % Calculate qtilde
        dat.qcum = cumsum(dat.dzqn_ref.*dat.dz);
        dat.qtilde = 1 - dat.qcum/dat.qcum(end);
        iz = find(dat.qtilde == 0, 1, 'first');
        dat.qtilde((iz+1):end) = nan;
        dat2.qcum = cumsum(dat2.dzqn_ref.*dat.dz);
        dat2.qtilde = 1 - dat2.qcum/dat2.qcum(end);
        iz = find(dat2.qtilde == 0, 1, 'first');
        dat2.qtilde((iz+1):end) = nan;
        
        % Calculate normalized mass fluxes
        dat.mass_flux_norm = ...
            dat.mass_flux_cond9999/max(dat.mass_flux_cond9999);
        dat2.mass_flux_norm = ...
            dat2.mass_flux_cond9999/max(dat2.mass_flux_cond9999);
        
        % Interpolate onto global grids
        dat.mass_flux_norm_p = ...
            interp1(dat.p, dat.mass_flux_norm, pgrid, 'linear', 0);
        dat.mass_flux_norm_q = ...
            interp1(dat.qtilde(~isnan(dat.qtilde)), ...
            dat.mass_flux_norm(~isnan(dat.qtilde)), qgrid, 'linear', 0);
        dat2.mass_flux_norm_p = ...
            interp1(dat2.p, dat2.mass_flux_norm, pgrid, 'linear', 0);
        dat2.mass_flux_norm_q = ...
            interp1(dat2.qtilde(~isnan(dat2.qtilde)), ...
            dat2.mass_flux_norm(~isnan(dat2.qtilde)), qgrid, 'linear', 0);
        
        % Calculate edge weights
        edges_p(ii,jj) = -trapz(pgrid, ...
            abs(dat.mass_flux_norm_p - dat2.mass_flux_norm_p)) ./ sqrt(...
            trapz(pgrid, dat.mass_flux_norm_p) * ...
            trapz(pgrid, dat2.mass_flux_norm_p));
        edges_q(ii,jj) = -trapz(qgrid, ...
            abs(dat.mass_flux_norm_q - dat2.mass_flux_norm_q)) ./ sqrt(...
            trapz(qgrid, dat.mass_flux_norm_q) * ...
            trapz(qgrid, dat2.mass_flux_norm_q));
        
        ii = ii+1;
        
    end
    
    jj = jj+1;
    
end

fprintf("p: %.2f (max %.2f, min %.2f)\n", ...
    mean(mean(edges_p(~isinf(edges_p)))), ...
    max(edges_p(~isinf(edges_p))), ...
    min(edges_p(~isinf(edges_p))));

fprintf("q: %.2f (max %.2f, min %.2f)\n", ...
    mean(mean(edges_q(~isinf(edges_q)))), ...
    max(edges_q(~isinf(edges_q))), ...
    min(edges_q(~isinf(edges_q))));
