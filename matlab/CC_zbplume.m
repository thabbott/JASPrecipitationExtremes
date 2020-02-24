%% CC_zbplume
% Calculate temperature profiles using a zero-buoyancy plume model
%
% Tristan Abbott // Massachusetts Institute of Technology // 07/03/2019
%
%%% Syntax
%   [T, qsat, p] = CC_zbplume(T0, qv0, p, refine, RH, entrain)
%
%%% Description
% Calculates temperature profiles set by a zero-buoyancy plume.
% The profiles start at a base state of $T = T0$, $q_v = qv0$, and $p =
% p(:,:,1)$ and are calculated by stepping upward 
% over the provided pressure grid to $p = p(:,:,end)$.
%
%%% Input Arguments
% *T - base temperature:*
% Temperatures at the base of each column. Must be in units of K, and my be
% provided as either a single value (for a single-column calculation) or as a 2D
% matrix (to simultaneously compute moist adiabats over multiple columns).
%
% *qv0 - base water vapor specific humidity:*
% Water vapor specific humidities at the base of each column. Must be in units
% of kg/kg and must match the shape of T. If this argument is larger than the
% saturation specific humidity at the base of a column, it will be set to the
% saturation specific humidity.
% 
% *p - pressure grid:*
% Pressure grid on which the moist adiabat is calculated, in units of Pa. The
% first two dimensions of this argument must match the dimensions of T. 
% This argument determines both the pressure range over which moist adiabatic 
% profiles are
% computed-- they begin in the state prescribed by T0, qv0, and p(:,:,1), and
% they end at the pressures given by p(:,:,end)-- and the size of the steps used
% to move upward along a moist adiabat. Increasing the resolution of this grid
% will therefore increase the resolution of the output profiles and decrease
% numerical errors, albeit at the expense of computation time.
%
% *refine - pressure grid refinement factor:*
% Refine the pressure grid before stepping upwards. If this factor is set to 1,
% successive steps move between successive levels of p; if this factor is larger
% than 1, that number of steps are taken between successive levels. Increasing
% the value of this argument will decrease numerical errors and increase the
% resolution of output profiles at the expense of increased computation time.
%
% *RH - environmental relative humidity:*
% Relative humidity of the environment surrounding the plume. Must be
% scalar.
%
% *entrain - plume entrainment rate:*
% Entrainment rate (1/m) of the zero-buoyancy plume. Must be scalar.
%
%%% Output Arguments
% *T - temperatures:*
% Temperature profiles, in units of K.
%
% *h - moist static energy:*
% Moist static energy profiles, in units of J/kg
%
% *qsat - saturation specific humidity:*
% Saturation specific humidity profiles.
%
% *p - pressure grid:*
% Pressure grid for the output profiles, in units of Pa.
%
% *z - height grid:*
% Height grid for the output profiles, in units of m.
%
%%% Source code
function [T, h, qstar, p, z] = CC_zbplume(T0, qv0, p, refine,...
    RH, entrain)  

    % Define constants
    global SAM_c_p;
    global SAM_L_c;
    global SAM_g;
    global SAM_R;
    
    % Set up pressure grid
    if refine > 1
        p = permute(p, [3 2 1]);
        p = interp1(p, (refine:refine*size(p,1))/refine);
        p = permute(p, [3 2 1]);
    end

    % Define function for saturation MSE
    sat_mse = @(T, z, p) SAM_c_p * T + SAM_g * z + SAM_L_c * ...
        SOM_qsat(T, p);
    
    % Create saturation specific humidity output
    qstar = zeros(size(p));
    qstar(:,:,1) = SOM_qsat(T0, p(:,:,1));

    % Reset supersaturated initial conditions to saturation
    m = logical(qv0 >= qstar(:,:,1));
    tmp1 = qstar(:,:,1);
    qv0(m) = tmp1(m);
    
    % Create height output
    z = zeros(size(p));
    z(:,:,1) = 0;
    
    % Calculate MSE at first level
    h = zeros(size(p));
    h(:,:,1) = SAM_c_p * T0 + SAM_g * z(:,:,1) + SAM_L_c * qv0;
    
    % Create temperature output
    T = zeros(size(p));
    T(:,:,1) = T0;
    
    % Step upwards over the grid
    fprintf('%06.2f/100', 0);
    for ii = 2:size(p,3)
        fprintf('\b\b\b\b\b\b\b\b\b\b%06.2f/100', ...
            100*ii/size(p,3));
        
        % Calculate height increment
        rho = p(:,:,ii-1) ./ (SAM_R * T(:,:,ii-1));
        dp = p(:,:,ii) - p(:,:,ii-1);
        dz = -dp ./ (rho * SAM_g);
        z(:,:,ii) = z(:,:,ii-1) + dz;
        
        % Calculate entrainment rate (0 below cloud base)
        entrain_act = entrain * ones(size(z(:,:,ii-1)));
        entrain_act(qv0 < qstar(:,:,ii-1)) = 0;
        
        % Calculate MSE increment
        dh = -entrain_act .* SAM_L_c .* (1 - RH) .* qstar(:,:,ii-1) .* dz;
        h(:,:,ii) = h(:,:,ii-1) + dh;
        
        % Invert for temperature
        for kk = 1:size(p,2)
            for jj = 1:size(p,1)
                T(jj,kk,ii) = fzero(...
                    @(tabs) h(jj,kk,ii) - ...
                    sat_mse(tabs, z(jj,kk,ii), p(jj,kk,ii)), ...
                    T(jj,kk,ii-1));
            end
        end
        
        % Calculate provisional saturation specific humidity
        qstar(:,:,ii) = SOM_qsat(T(:,:,ii), p(:,:,ii));
        
        % Identify unsaturated parcels and correct temperature
        mask = qstar(:,:,ii) > qv0;
        if sum(mask) > 0
            T_all = T(:,:,ii);
            T_unsat = (h(:,:,ii) - SAM_g * z(:,:,ii) - SAM_L_c * qv0) / ...
                SAM_c_p;
            T_all(mask) = T_unsat(mask);
            T(:,:,ii) = T_all;
        end
        
        % Calculate new saturation specific humidity
        qstar(:,:,ii) = SOM_qsat(T(:,:,ii), p(:,:,ii));

    end
    fprintf('\n');

end
