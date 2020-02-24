%% CC_pseudoAdiabat
% Calculate temperature and water vapor profiles along a moist pseudo-adiabat.
%
% Tristan Abbott // Massachusetts Institute of Technology // 01/16/2016
%
%%% Syntax
%   [T, qv, dpqv, qsat, p] = CC_pseudoAdiabat(T0, qv0, p, refine)
%
%%% Description
% Calculates temperature and water vapor profiles along a moist pseudo-adiabat.
% The pseudo-adiabats start at a base state of $T = T0$, $q_v = qv0$, and $p =
% p(:,:,1)$ and are calculated by stepping upward 
% over the provided pressure grid to $p = p(:,:,end)$.
%
% Below the saturation specific humidity, the profile follows a dry adiabat:
%
% $$
% T(p) = T_0 \left (\frac{p}{p_0} \right)^{R/c_p}
% $$
% 
% where $R = (1 - q_v)R_a + q_v R_v$ and $c_p = (1 - q_v)c_{pa} + q_v c_{pv}$.
% Over this range, $q_v = q_{v0}$ and $\partial_p q_v = 0$.
%
% Once the profile reaches a temperature and pressure where $q^* = q_v$, the
% profile follows a moist pseudo-adiabat. The temperature profile is built up
% using the equation for $dT/dp$ given in Pierrehumbert (2010), 
% _Principles of Planetary Climate_, equation (2.33), assuming the dilute limit
% (i.e. $r^* \approx q^*$). The moisture profile is build up based on the
% Clausius-Clapeyron relation for an ideal gas:
%
% $$
% \frac{de^*}{dT} = \frac{Le^*}{R_v T^2}.
% $$
%
% By assuming the dilute limit, where $q^* \approx e^* R_a R_v^{-1} p^{-1}$, the
% Clausius-Clapeyron relation can be rewritten in terms of $q^*p$ and solved for
% a pressure derivative of the saturation specific humidity:
%
% $$
% \frac{dq^*}{dp} = \frac{q^*}{p} \left [ \frac{dT}{dp} \frac{L q^* p}{R_v T^2}
% - 1 \right ].
% $$
%
% This function assumes that all condensate falls out as water, and $L$
% therefore represents the latent heat of condensation.
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
%%% Output Arguments
% *T - temperatures:*
% Temperatures along the moist adiabat, in units of K.
%
% *qv - water vapor:*
% Water vapor specific humidities along the moist adiabat, in units of kg/kg.
%
% *dpqv - water vapor lapse rates:*
% Pressure derivative of water vapor specific humidities along the moist
% adiabat, in units of kg/kg/Pa.
%
% *qsat - saturation specific humidity:*
% Saturation specific humidity at each level along the moist adiabat.
%
% *p - pressure grid:*
% Pressure grid for the output profiles, in units of Pa.
%
%%% <../test/html/CC_pseudoAdiabat_test.html Tests>
%
%%% Source code
function [T, qv, dpqv, qsat, p] = CC_pseudoAdiabat(T0, qv0, p, refine)  

    % Define constants
    global SAM_R;
    global SAM_R_v;
    global SAM_c_p;
    global SAM_c_pv;
    global SAM_L_c;
    
    % Calculate dry adiabatic R/cp
    Rcp = ((1 - qv0)*SAM_R + qv0*SAM_R_v) ./...
        ((1 - qv0)*SAM_c_p + qv0*SAM_c_pv);

    % Define anonymous functions for moist temperature derivative...
    dTdp = @(T, p, qv) ((T .* SAM_R) ./ (p * SAM_c_p)) .*...
        ((1 + (SAM_L_c.*qv)./(SAM_R.*T)) ./ ...
        (1 + qv.*(SAM_c_pv./SAM_c_p + (SAM_L_c./(SAM_c_p.*T)).*...
        ((SAM_L_c./(SAM_R_v.*T)) - 1))));
    % ... dry temperature derivative ...
    dTdpd = @(T, p) Rcp.*T./p; 
    % ... and water vapor derivative...
    dqvdp = @(T, p, qv, dTdp) (qv./p) .* ...
        ((dTdp .* SAM_L_c .* p) ./ (SAM_R_v .* T.^2) - 1);

    % Set up pressure grid
    if refine > 1
        p = permute(p, [3 2 1]);
        p = interp1(p, (refine:refine*size(p,1))/refine);
        p = permute(p, [3 2 1]);
    end
    dp = diff(p, 1, 3);

    % Check for base water vapor values that are above saturation
    qsat = zeros(size(p));
    qsat(:,:,1) = SAM_qsatWater(T0, p(:,:,1));
    m = logical(qv0 >= qsat(:,:,1));
    tmp1 = qsat(:,:,1);
    qv0(m) = tmp1(m);
    
    % Calculate saturation vapor pressure at next level
    dpqv = zeros(size(p));
    dpT = zeros(size(T0));
    tmp1 = dTdp(T0, p(:,:,1), qv0);
    dpT(m) = tmp1(m);
    tmp1 = dTdpd(T0, p(:,:,1));
    dpT(~m) = tmp1(~m);
    dpqv(:,:,1) = dqvdp(T0, p(:,:,1), qsat(:,:,1), dpT);
    qsat(:,:,2) = qsat(:,:,1) + dpqv(:,:,1).*dp(:,:,1);
    
    % Create output arrays and fill at first pressure level
    % Temperature
    T = zeros(size(p));
    T(:,:,1) = T0;
    % Water vapor
    qv = zeros(size(p));
    qv(:,:,1) = qv0;
    % Lapse rates
    tmp1 = dqvdp(T0, p(:,:,1), qv0, dpT);
    tmp1(~m) = 0;
    dpqv(:,:,1) = tmp1;
    
    % Step upwards over the grid
    for ii = 2:size(p,3)

        % Calculate temperature at the current step
        T(:,:,ii) = T(:,:,ii-1) + dpT.*dp(:,:,ii-1);
        
        % Calculate water vapor at current step
        qv(:,:,ii) = qv(:,:,ii-1) + dpqv(:,:,ii-1).*dp(:,:,ii-1);

        % Find points that have transitioned to moist adiabats
        m = logical(qv(:,:,ii) >= qsat(:,:,ii));
        tmp1 = qsat(:,:,ii);
        tmp2 = qv(:,:,ii);
        tmp2(m) = tmp1(m);
        qv(:,:,ii) = tmp2;

        % Calculate lapse rates at current step
        tmp1 = dTdp(T(:,:,ii), p(:,:,ii), qv(:,:,ii));
        dpT(m) = tmp1(m);
        tmp1 = dTdpd(T(:,:,ii), p(:,:,ii));
        dpT(~m) = tmp1(~m);
        
        % Calculate saturation specific humidity at next step
        if ii < size(p,3)
            dpqv(:,:,ii) = dqvdp(T(:,:,ii), p(:,:,ii), qsat(:,:,ii), dpT);
            qsat(:,:,ii+1) = qsat(:,:,ii) + dpqv(:,:,ii).*dp(:,:,ii);
        end

        % Calculate water vapor lapse rates at current step
        tmp1 = dqvdp(T(:,:,ii), p(:,:,ii), qv(:,:,ii), dpT);
        tmp1(~m) = 0;
        dpqv(:,:,ii) = tmp1;

    end

end
