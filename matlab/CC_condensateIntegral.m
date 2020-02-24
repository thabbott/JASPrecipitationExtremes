%% CC_condensateIntegral
% Calculate a condensate integral using a moist pseudo-adiabat based on the
% Clausius-Clapeyron relation.
%
% Tristan Abbott // Massachusetts Institute of Technology // 02/01/2016
%
%%% Syntax
%  [C, dpqv, p] = CC_condensateIntegral(T0, qv0, w, p, refine)
%  [...] = CC_condensateIntegral(..., 'domain', 'updrafts')
%  [...] = CC_condensateIntegral(..., 'mask', m)
%
%%% Description
%
% Calculate a condensate integral over a moist pseudo-adiabat defined based on
% the Clausius-Clapeyron relation for ideal gases. The definition of the
% pseudoadiabat is given in <CC_pseudoAdiabat.html CC_pseudoAdiabat>. The
% condensate integral is
%
% $$
% C = \int_{p_0}^{p_1} \frac{\partial q_v}{\partial p}
% \omega(p) g^{-1} \mathrm{d}\,p.
% $$
%
% $q_v$ is the saturation specific humidity, and its pressure derivative is taken
% along a moist pseudo-adiabat. $\omega$ is the vertical velocity in pressure
% coordinates.
%
% This integral is the same as the one computed in
% <SOM_condensateIntegral.html
% SOM_condensateIntegral>; here, hydrostatic balance is used to express it in
% pressure rather than height coordinates.
%
%
%%% Input Arguments
% *T0 - base temperature:*
% Temperatures at the base of each column; that is, at the pressures given in
% p(:,:,1). Must be in units of K, and my be
% provided as either a single value (for a single-column calculation) or as a 2D
% matrix (to simultaneously compute moist adiabats over multiple columns).
%
% *qv0 - base water vapor specific humidity:*
% Water vapor specific humidities at the base of each column; that is, at the
% pressures given in p(:,:,1). Must be in units
% of kg/kg and must match the shape of T. If this argument is larger than the
% saturation specific humidity at the base of a column, it will be set to the
% saturation specific humidity.
%
% *w - pressure velocity/mass flux profiles:*
% Mass flux profiles in each column, in units of Pa/s. The first two 
% dimensions of this argument must match the dimensions of T.
%
% *p - pressure grid:*
% Pressure grid on which the integral is evaluated, in units of Pa. The
% first two dimensions of this argument must match the dimensions of T, and the
% third dimension must match the third dimension of omega. This argument
% determines the pressure range over which moist adiabatic profiles are
% computed: they begin in the state prescribed by T0, qv0, and p(:,:,1), and
% they end at the pressures given by p(:,:,end).
%
% *refine - refinement factor:*
% Refine the pressure grid before stepping upwards. If this factor is set to 1,
% successive steps move between successive levels of p; if this factor is larger
% than 1, that number of steps are taken between successive levels. Increasing
% the value of this argument will decrease numerical errors and increase the
% resolution of output profiles at the expense of increased computation time.
% The condensate integral is computed using the refined grid, and omega is
% interpolated linearly onto it before doing so.
%
%%% Parameters
% *domain - domain criterion for the integral:*
% Specify a criterion that sections of the profiles must meet in order to
% be included in the integral. Options are 'all' (include the entire
% profile) and 'updrafts' (include only sections of the profile where the
% vertical velocity is upwards).
%
% *mask - integral domain mask:* a mask (logical array) that can exclude 
% elements from condensation integrals. Allowable shapes are the same as 
% for z, and elements of integrands that are not masked will be set to 0
% before integrating.
%
%%% Output Arguments
% *C - condensate integral:*
% Value of the condensate integral over each column, in units of kg of water
% per square meter per second.
%
% *dpqv - pressure derivative of water vapor:*
% Pressure derivative of water vapor at each point on the refined grid.
% Units are kg/kg/Pa.
%
% *omega - pressure velocity:*
% Pressure velocity at each point on the refined grid. Units are Pa/s.
%
% *p - refined pressure grid (Pa).*
%
%%% <../test/html/CC_condensateIntegral_test.html Tests>
%
%%% Source code
function [C, dpqv, p] = CC_condensateIntegral(T0, qv0, w, p, refine, ...
    varargin)  
	
    global SAM_g;
    
    % Check for parameters
    domain = 0;
    masked = 0;
    mask = [];
    if numel(varargin) ~= 0
        for i = 2:2:numel(varargin)
            switch lower(varargin{i-1})
                case 'domain'
                    switch lower(varargin{i})
                        case 'updrafts'
                            domain = 1;
                        case 'all'
                        otherwise
                            error(['Unrecognized value for parameter',...
                                ' "domain": %s'], varargin{i});
                    end
                case 'mask'
                    masked = 1;
                    mask = varargin{i};
                otherwise
                    error('Unrecognized parameter: %s', varargin{i-1});
            end
        end
    end

    % Set up the pressure grid
    if refine > 1
        p = permute(p, [3 2 1]);
        p = interp1(p, (refine:refine*size(p,1))/refine);
        p = permute(p, [3 2 1]);
    end

    % Calculate moist adiabatic profiles
    [~,~,dpqv,~,p] = CC_pseudoAdiabat(T0, qv0, p, 1);
    
    % Interpolate w onto the refined grid
    w = permute(w, [3 2 1]);
    w = interp1(w, (refine:refine*size(w,1))/refine);
    w = permute(w, [3 2 1]);
    
    % Multiply moisture derivative by pressure velocity
    dpqvw = dpqv.*w;
    
    % Zero part of the integral if needed
    if domain == 1
        dpqvw(logical(w > 0)) = 0;
    end
    if masked == 1
        dpqvw(~mask) = 0; 
    end
    
    % Compute integral with trapezoidal approximation along z
    integrand = (p(:,:,2:end) - p(:,:,1:end-1)) .* ...
        (dpqvw(:,:,2:end) + dpqvw(:,:,1:end-1));
    C = sum(integrand, 3) ./ (2*SAM_g);

end
