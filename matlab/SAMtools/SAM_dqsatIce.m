%% SAM_dqsatIce
% Temperature derivative of saturation specific humidity over ice
%
%%% Syntax
%   dq = SAM_dqsatIce(T)
%
%%% Description
% Calculates the temperature derivative of saturation specific humidity over 
% ice with a
% polynomial fit. The polynomial is based on Flatau, Walko, and
% Cotton (1992): "Polynomial Fits to Saturation Vapor Pressure" and is
% identical to the one used in the System for Atmospheric Modeling, version
% 6.10.8.
%
%%% Input Arguments
% *T - temperature (K):*
% May be either scalar or non-scalar. If non-scalar, the output has the same
% size and shape as the input.
%
% *p - pressure (Pa):*
% Ambient pressure. Must have the same dimensions as T.
%
%%% Output Arguments
% *dq - temperature derivative of saturation specific humidity (nondim):*
% Temperature derivative of saturation specific humidity, 
% in units of kg water per total kg per Kelvin.
%
%%% <../test/html/SAM_dqsatIce_test.html Tests>

function dq = SAM_dqsatIce(T, p)  

    dq = 0.622 * SAM_dpsatIce(T) ./ p;

end
