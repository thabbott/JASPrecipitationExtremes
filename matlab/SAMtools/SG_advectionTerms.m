%% SG_advectionTerms
% Approximate the advection terms for a field located at scalar
% cell centers.
%
% Tristan Abbott // Massachusetts Institute of Technology // 11/14/2017
%
%%% Syntax
%   advw = SG_advectionTerms(grid, 's', 'u', 'v', 'w');
%
%%% Description
% Approximates the advection terms-- u ds/dx and so on for a field s-- by 
% using a centered first-order finite difference to compute ds/dz at
% vertical velocity levels. This type of discretization generally does not
% have good conservation properties, and SG_advectScalar should be used if
% the goal is to compute total advective tendencies. One additional
% important difference between this function and SG_advectScalar is the
% sign of the result-- this function returns the value of the advection
% terms (e.g. u ds/dx) whereas SG_advectScalar returns advective
% tendencies, which have the opposite sign.
%
%%% Input Arguments
% *grid - SAM C-grid struct:*
% A struct containing grid information generated by SG_grid and fields
% added by SG_addVar.
%
% *s, u, v, w - names of scalar and velocity:*
% String keys that can be used to look up scalar and velocity fields held
% in the grid input struct.
%
%%% Output Arguments
%
% *ax, ay, az - advection terms:*
% Arrays containing the values of the advection terms at velocity points.
%

function [ax, ay, az] = SG_advectionTerms(grid, s, u, v, w)

    % U direction
    ds = zeros(size(grid.(u)));
    dx = grid.x(2) - grid.x(1);
    ds(:,:,1) = (grid.(s)(:,:,1) - grid.s(:,:,end))/dx;
    ds(:,:,end) = ds(:,:,1);
    ds(:,:,2:end-1) = (grid.(s)(:,:,2:end) - grid.(s)(:,:,1:end-1))/dx;
    ax = grid.(u).*ds;
    
    % V direction
    ds = zeros(size(grid.(v)));
    dy = grid.y(2) - grid.y(1);
    ds(:,1,:) = (grid.(s)(:,1,:) - grid.s(:,end,:))/dy;
    ds(:,end,:) = ds(:,1,:);
    ds(:,2:end-1,:) = (grid.(s)(:,2:end,:) - grid.(s)(:,1:end-1,:))/dy;
    ay = grid.(v).*ds;
    
    % Vertical
    ds = zeros(size(grid.(w)));
    ds(2:end-1,:,:) = (grid.(s)(2:end,:,:) - grid.(s)(1:end-1,:,:)) ./...
        (grid.z(2:end) - grid.z(1:end-1));
    az = grid.(w).*ds;

end