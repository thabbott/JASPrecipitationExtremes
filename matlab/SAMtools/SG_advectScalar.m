%% SG_grid
% Compute advective tendencies for a value held at scalar points
% on the SAM C-grid.
%
% Tristan Abbott // Massachusetts Institute of Technology // 11/14/2017
%
%%% Syntax
%   [ax, ay, az] = SG_advectScalar(grid, 's', 'u', 'v', 'w', 'upwind1');
%
%%% Description
% Computes tendencies from advective fluxes in x, y, and z directions and
% returns each tendency separately.
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
% *advopt - advection scheme:*
% String containing the advection scheme used to compute tendencies.
% Supported options are 'upwind1' (first-order upwind scheme),
% 'adv_mpdata' (SAM advection option)
%
%%% Output Arguments
%
% *ax, ay, az - advective tendencies:*
% Arrays containing advective tendencies from flow in the x, y, and z
% directions. Note that these are advective tendencies from the flux form
% of the advection equation, so ax = -1/rho * d/dx (rho u s), and the
% vertical advection term w ds/dz is contained in ax + ay.
%
function [ax, ay, az] = SG_advectScalar(grid, s, u, v, w, advopt)
    
    switch lower(advopt)
        case 'upwind1'
            [ax, ay, az] = upwind1(grid, s, u, v, w);
        case 'adv_mpdata'
            [ax, ay, az] = adv_mpdata(grid, s, u, v, w);
        otherwise
            error('SG_advectScalar: unknown advection scheme %s',...
                lower(advopt))
    end
    
end

function [ax, ay, az] = upwind1(grid, s, u, v, w)
    uuu = zeros(size(grid.(u)));
    vvv = zeros(size(grid.(v)));
    www = zeros(size(grid.(w)));
    nx = size(grid.(s), 3);
    ny = size(grid.(s), 2);
    nz = size(grid.(s), 1);  
        
    % Compute metric scaling factor
    dz = grid.zi(2:end) - grid.zi(1:end-1);
    dx = grid.x(2) - grid.x(1);
    dy = grid.y(2) - grid.y(1);
    
    % Create 3d density grids
    rho3d = repmat(grid.rho, [1 ny nx]);
    rhoi3d = repmat(grid.rhoi, [1 ny nx]);
    
    % Upwind interpolation for u*s
    ii = 1;
    uuu(:,:,ii) = ...
        max(0, grid.(u)(:,:,ii)).*grid.(s)(:,:,end) + ...
        min(0, grid.(u)(:,:,ii)).*grid.(s)(:,:,ii);
    uuu(:,:,end) = uuu(:,:,ii);
    for ii = 2:nx
        uuu(:,:,ii) = ...
            max(0, grid.(u)(:,:,ii)).*grid.(s)(:,:,ii-1) + ...
            min(0, grid.(u)(:,:,ii)).*grid.(s)(:,:,ii);
    end
    
    % Upwind interpolation for v*s
    jj = 1;
    vvv(:,jj,:) = ...
        max(0, grid.(v)(:,jj,:)).*grid.(s)(:,end,:) + ...
        min(0, grid.(v)(:,jj,:)).*grid.(s)(:,jj,:);
    vvv(:,end,:) = vvv(:,jj,:);
    for jj = 2:ny
        vvv(:,jj,:) = ...
        max(0, grid.(v)(:,jj,:)).*grid.(s)(:,jj-1,:) + ...
        min(0, grid.(v)(:,jj,:)).*grid.(s)(:,jj,:);
    end
    
    % Upwind interpolation for w*s
    for kk = 2:nz
        www(kk,:,:) = ...
            max(0, grid.(w)(kk,:,:)).*grid.(s)(kk-1,:,:) + ...
            min(0, grid.(w)(kk,:,:)).*grid.(s)(kk,:,:);
    end
    
    % Compute tendency from U convergence
    ax = 1/dx * (uuu(:,:,1:end-1) - uuu(:,:,2:end));
    
    % Compute tendency from V convergence
    ay = 1/dy * (vvv(:,1:end-1,:) - vvv(:,2:end,:));
    
    % Compute tendency from W convergence
    idz = 1./(repmat(dz, [1 ny nx]));
    az = idz .* (rhoi3d(1:end-1,:,:).*www(1:end-1,:,:) - ...
        rhoi3d(2:end,:,:).*www(2:end,:,:)) ./ rho3d;
    
end

function [ax, ay, az] = adv_mpdata(grid, ss, us, vs, ws)
    
    % Compute dimensions
    d = size(grid.(ss));
    nx = d(3); ny = d(2); nzm = d(1); nz = nzm+1;

    % Allocate arrays
    f = zeros(nx+6,ny+6,nzm);
    u = zeros(nx+5,ny+4,nzm);
    v = zeros(nx+4,ny+5,nzm);
    w = zeros(nx+4,ny+4,nz);
    mx = zeros(nx+2,ny+2,nzm);
    mn = zeros(nx+2,ny+2,nzm);
    uuu = zeros(nx+5,ny+4,nzm);
    vvv = zeros(nx+4,ny+5,nzm);
    www = zeros(nx+4,ny+4,nz);
    iadz = 1./grid.adz;
    irho = 1./grid.rho;
    irhow = 1./grid.rhoi;

    % Exchange values on input boundaries
    f(4:nx+3,4:ny+3,:) = permute(grid.(ss), [3 2 1]);
    u(3:nx+2,3:ny+2,:) = permute(grid.(us)(:,:,1:nx), [3 2 1]);
    v(3:nx+2,3:ny+2,:) = permute(grid.(vs)(:,1:ny,:), [3 2 1]);
    w(3:nx+2,3:ny+2,:) = permute(grid.(ws), [3 2 1]);
    
    % General case
    if (nx > 2 && ny > 2)
        f(1:3,:,:) = f(nx+1:nx+3,:,:);
        f(nx+4:end,:,:) = f(4:6,:,:);
        f(:,1:3,:) = f(:,ny+1:ny+3,:);
        f(:,ny+4:end,:) = f(:,4:6,:);
        
        u(1:2,:,:) = u(nx+1:nx+2,:,:);
        u(nx+3:end,:,:) = u(3:5,:,:);
        u(:,1:2,:) = u(:,ny+1:ny+2,:);
        u(:,ny+3:end,:) = u(:,3:4,:);
        
        v(1:2,:,:) = v(nx+1:nx+2,:,:);
        v(nx+3:end,:,:) = v(3:4,:,:);
        v(:,1:2,:) = v(:,ny+1:ny+2,:);
        v(:,ny+3:end,:) = v(:,3:5,:);
        
        w(1:2,:,:) = w(nx+1:nx+2,:,:);
        w(nx+3:end,:,:) = w(3:4,:,:);
        w(:,1:2,:) = w(:,ny+1:ny+2,:);
        w(:,ny+3:end,:) = w(:,3:4,:);
    % Special case: column
    elseif (nx == 1 && ny == 1)
        f = repmat(f(4,4,:), [size(f,1) size(f,2) 1]);
        u = repmat(u(3,3,:), [size(u,1) size(u,2) 1]);
        v = repmat(v(3,3,:), [size(v,1) size(v,2) 1]);
        w = repmat(w(3,3,:), [size(w,1) size(w,2) 1]);
    % Otherwise, can't do it (at least for now)
    else
        error('SG_advectScalar:adv_mpdata:dimensions', ...
            ['Strange number of dimensions in input data: ',...
            'nx = ', str(nx), ', ny = ', str(ny), '. More than ',...
            'a single column but not enough for horizontal advection']);
    end

    % Define helper for reshaping reference profiles
    mangle = @(p,i,j,k) repmat(reshape(p(k), [1 1 numel(k)]),...
                               [numel(i) numel(j) 1]);

    % Define inline functions
    andiff = @(x1,x2,a,b) (abs(a)-a.*a.*b).*0.5.*(x2-x1);
    across = @(x1,a1,a2) 0.03125.*a1.*a2.*x1;
    pp = @(y) max(0,y);
    pn = @(y) -min(0,y);

    % Other parameters
    eps = 1e-10;

    % Debugging
    ip = 91; jp = 14; kp = 27;
    doprint = 0;

    % Skipping if (dowallx), if (dowally)

    if doprint == 1
        printf(f,ip,jp,kp);
        ind = find(f > 318.95983 & f < 318.95984);
        [i,j,k] = ind2sub(size(f), ind);
        printf(f,i-3,j-3,k);
    end

    % if (nonos) then
    % Indexing base f, shift (-2,-2,0) for mx, mn
    % for kk = 1:nzm
    kk = 1:nzm;
        k = kk;
        kc = min(nzm,k+1);
        kb = max(1, k-1);
        % for jj = 0:ny+1
        jj = 0:ny+1;
            j = jj+3;
            jb = j-1;
            jc = j+1;
            % for ii = 0:nx+1
            ii = 0:nx+1;
                i = ii+3;
                ib = i-1;
                ic = i+1;
                mx(i-2,j-2,k)=max(f(ib,j,k),max(f(ic,j,k),max(f(i,jb,k),...
                              max(f(i,jc,k),max(f(i,j,kb),max(f(i,j,kc),...
                              f(i,j,k)))))));
                mn(i-2,j-2,k)=min(f(ib,j,k),min(f(ic,j,k),min(f(i,jb,k),...
                              min(f(i,jc,k),min(f(i,j,kb),min(f(i,j,kc),...
                              f(i,j,k)))))));
            % end % for ii = 0:nx+1
        % end % for jj = 0:ny+1
    % end % for kk = 1:nzm
    % endif (nonos)

    if doprint == 1
        printmx(mx,ip,jp,kp);
        printmn(mn,ip,jp,kp);
    end

    % Indexing base u, shift (+1,+1,0) for f
    % for kk = 1:nzm
    kk = 1:nzm;
        k = kk;
        % for jj = -1:ny+2
        jj = -1:ny+2;
            j = jj+2;
            % for ii = -1:nx+3
            ii = -1:nx+3;
                i = ii+2;
                uuu(i,j,k) = max(0,u(i,j,k)).*f(i-1+1,j+1,k) + ...
                             min(0,u(i,j,k)).*f(i+1,j+1,k);
            % end % for ii = -1:nx+3
        % end % for jj = -1:ny+2
    % end % for kk = 1:nzm
    if doprint == 1
        printuuu(uuu,ip,jp,kp);
    end
    
    % Indexing base v, shift (+1,+1,0) for f
    % for kk = 1:nzm
    kk = 1:nzm;
        k = kk;
        % for jj = -1:ny+3
        jj = -1:ny+3;
            j = jj+2;
            % for ii = -1:nx+2
            ii = -1:nx+2;
                i = ii+2;
                vvv(i,j,k) = max(0,v(i,j,k)).*f(i+1,j-1+1,k) + ...
                             min(0,v(i,j,k)).*f(i+1,j+1,k);
            % end % for ii = -1:nx+3
        % end % for jj = -1:ny+2
    % end % for kk = 1:nzm
    if doprint == 1
        printvvv(vvv,ip,jp,kp);
    end
     
    % Indexing base w, shift (+1,+1,0) for f
    % for kk = 1:nzm
    kk = 1:nzm;
        k = kk;
        kb = max(1,k-1);
        % for jj = -1:ny+2
        jj = -1:ny+2;
            j = jj+2;
            % for ii = -1:nx+2
            ii = -1:nx+2;
                i = ii+2;
                www(i,j,k) = max(0,w(i,j,k)).*f(i+1,j+1,kb) + ...
                             min(0,w(i,j,k)).*f(i+1,j+1,k);
            % end % for ii = -1:nx+3
        % end % for jj = -1:ny+2
    % end % for kk = 1:nzm
    if doprint == 1
        printwww(www,ip,jp,kp);
    end

    %===
    % LINE 143
    %===
    % Indexing base uvw, shift (+1,+1,0) for f
    % for kk = 1:nzm
    kk = 1:nzm;
        k = kk;
        % for jj = -1:ny+2
        jj = -1:ny+2;
            j = jj+2;
            % for ii = -1:nx+2
            ii = -1:nx+2;
                i = ii+2;
                f(i+1,j+1,k) = f(i+1,j+1,k) - (uuu(i+1,j,k) - uuu(i,j,k) + ...
                                    vvv(i,j+1,k) - vvv(i,j,k) + ...
                                    (www(i,j,k+1) - www(i,j,k)).*...
                                    mangle(iadz,i,j,k)).* ...
                                    mangle(irho,i,j,k);
            % end % for ii = -1:nx+2
        % end % for jj = -1:ny+2
    % end % for kk = 1:nzm
    if doprint == 1
        printf(f,ip,jp,kp);
    end

    % Indexing base uvw, shift (+1,+1,0) for f
    % for kk = 1:nzm
    kk = 1:nzm;
        k = kk;
        kc = min(nzm,k+1);
        kb = max(1,k-1);
        % for jj = 0:ny+1
        jj = 0:ny+1;
            j = jj+2;
            jb = j-1;
            jc = j+1;
            % for ii = 0:nx+2
            ii = 0:nx+2;
                i = ii+2;
                ib = i-1;
        dd = 2./(mangle(kc,i,j,k)-mangle(kb,i,j,k))./mangle(grid.adz,i,j,k);
                uuu(i,j,k)=andiff(f(ib+1,j+1,k),f(i+1,j+1,k),u(i,j,k),...
                    mangle(irho,i,j,k)) ...
                    -(across(f(ib+1,jc+1,k)+f(i+1,jc+1,k)-f(ib+1,jb+1,k)-f(i+1,jb+1,k), ...
                    u(i,j,k), v(ib,j,k)+v(ib,jc,k)+v(i,jc,k)+v(i,j,k)) ...
                    +across(dd.*(f(ib+1,j+1,kc)+f(i+1,j+1,kc)-f(ib+1,j+1,kb)-f(i+1,j+1,kb)), ...
                    u(i,j,k), w(ib,j,k)+w(ib,j,kc)+w(i,j,k)+w(i,j,kc))) .*...
                    mangle(irho,i,j,k);
            % end % for ii = 0:nx+2
        % end % for jj = 0:ny+1
    % end % for kk = 1:nzm
    if doprint == 1
        printuuu(uuu,ip,jp,kp);
    end

    % Indexing base uvw, shift (+1,+1,0) for f
    % for kk = 1:nzm
    kk = 1:nzm;
        k = kk;
        kc = min(nzm,k+1);
        kb = max(1,k-1);
        % for jj = 0:ny+2
        jj = 0:ny+2;
            j = jj+2;
            jb = j-1;
            % for ii = 0:nx+1
            ii = 0:nx+1;
                i = ii+2;
                ib = i-1;
                ic = i+1;
        dd = 2./(mangle(kc,i,j,k)-mangle(kb,i,j,k))./mangle(grid.adz,i,j,k);
                vvv(i,j,k)=andiff(f(i+1,jb+1,k),f(i+1,j+1,k),v(i,j,k),...
                    mangle(irho,i,j,k)) ...
                    -(across(f(ic+1,jb+1,k)+f(ic+1,j+1,k)-f(ib+1,jb+1,k)-f(ib+1,j+1,k), ...
                    v(i,j,k), u(i,jb,k)+u(i,j,k)+u(ic,j,k)+u(ic,jb,k)) ...
                    +across(dd.*(f(i+1,jb+1,kc)+f(i+1,j+1,kc)-f(i+1,jb+1,kb)-f(i+1,j+1,kb)), ...
                    v(i,j,k), w(i,jb,k)+w(i,j,k)+w(i,j,kc)+w(i,jb,kc))) .*...
                    mangle(irho,i,j,k);
            % end % for ii = 0:nx+1
        % end % for jj = 0:ny+2
    % end % for kk = 1:nzm
    if doprint == 1
        printvvv(vvv,ip,jp,kp);
    end

    % Indexing base uvw, shift (+1,+1,0) for f
    % for kk = 1:nzm
    kk = 1:nzm;
        k = kk;
        kb = max(1,k-1);
        irhow(k) = 1./(grid.rhoi(k).*grid.adz(k));
        % for jj = 0:ny+1
        jj = 0:ny+1;
            j = jj+2;
            jb = j-1;
            jc = j+1;
            % for ii = 0:nx+1
            ii = 0:nx+1;
                i = ii+2;
                ib = i-1;
                ic = i+1;
                www(i,j,k)=andiff(f(i+1,j+1,kb),f(i+1,j+1,k),w(i,j,k),...
                    mangle(irhow,i,j,k)) ...
                   -(across(f(ic+1,j+1,kb)+f(ic+1,j+1,k)-f(ib+1,j+1,kb)-f(ib+1,j+1,k), ...
                   w(i,j,k), u(i,j,kb)+u(i,j,k)+u(ic,j,k)+u(ic,j,kb)) ...
                   +across(f(i+1,jc+1,k)+f(i+1,jc+1,kb)-f(i+1,jb+1,k)-f(i+1,jb+1,kb), ...
                   w(i,j,k), v(i,j,kb)+v(i,jc,kb)+v(i,jc,k)+v(i,j,k))) .* ...
                   mangle(irho,i,j,k);
            % end % for ii = 0:nx+1
        % end % for jj = 0:ny+1
    % end % for kk = 1:nzm
    if doprint == 1
        fprintf(1, 'irhow(%d) = %.6f\n', kp, irhow(kp));
        fprintf(1, 'irho(%d) = %.6f\n', kp, irho(kp));
        i=ip+2; j=jp+2; k=kp; kb=k-1; ib=i-1; ic=i+1; jb=j-1; jc=j+1;
        fprintf(1, 'andiff = %.9f\n', andiff(f(i+1,j+1,kb),f(i+1,j+1,k),w(i,j,k),irhow(k)));
        fprintf(1, 'across = %.9f\n', across(f(ic+1,j+1,kb)+f(ic+1,j+1,k)-f(ib+1,j+1,kb)-f(ib+1,j+1,k), w(i,j,k), u(i,j,kb)+u(i,j,k)+u(ic,j,k)+u(ic,j,kb)));
        fprintf(1, 'across = %.9f\n', across(f(i+1,jc+1,k)+f(i+1,jc+1,kb)-f(i+1,jb+1,k)-f(i+1,jb+1,kb),w(i,j,k), v(i,j,kb)+v(i,jc,kb)+v(i,jc,k)+v(i,j,k)));
        printwww(w,ip,jp,kp);
        printwww(www,ip,jp,kp);
    end
    
    www(:,:,1) = 0;

    %===
    % LINE 213
    %===
    % if (nonos) then
    % Indexing base f, shift (-2,-2,0) for mx,mn
    for kk = 1:nzm
        k = kk;
        kc = min([nzm,k+1]);
        kb = max([1,k-1]);
        for jj = 0:ny+1
            j = jj+3;
            jb = j-1;
            jc = j+1;
            for ii = 0:nx+1
                i = ii+3;
                ib = i-1;
                ic = i+1;
                mx(i-2,j-2,k)=max([f(ib,j,k),f(ic,j,k),f(i,jb,k), ...
                          f(i,jc,k),f(i,j,kb),f(i,j,kc),f(i,j,k),mx(i-2,j-2,k)]);
                mn(i-2,j-2,k)=min([f(ib,j,k),f(ic,j,k),f(i,jb,k), ...
                          f(i,jc,k),f(i,j,kb),f(i,j,kc),f(i,j,k),mn(i-2,j-2,k)]);
            end % for ii = 0:nx+1
        end % for jj = 0:ny+1
    end % for kk = 1:nzm
    if doprint == 1
        printmx(mx,ip,jp,kp);
        printmn(mn,ip,jp,kp);
    end

    % Indexing base u, shift (-1,-1,0) for mx,mn, (+1,+1,0) for f
    % for kk = 1:nzm
    kk = 1:nzm;
        k = kk;
        kc = min(nzm,k+1);
        % for jj = 0:ny+1
        jj = 0:ny+1;
            j = jj+2;
            jc = j+1;
            % for ii = 0:nx+1
            ii = 0:nx+1;
                i = ii+2;
                ic = i+1;
                mx(i-1,j-1,k)=mangle(grid.rho,i,j,k).*...
                          (mx(i-1,j-1,k)-f(i+1,j+1,k))./ ...
                          (pn(uuu(ic,j,k)) + pp(uuu(i,j,k))+ ...
                          pn(vvv(i,jc,k)) + pp(vvv(i,j,k))+ ...
                          mangle(iadz,i,j,k).*...
                          (pn(www(i,j,kc)) + pp(www(i,j,k)))+eps);  
                mn(i-1,j-1,k)=mangle(grid.rho,i,j,k).*...
                          (f(i+1,j+1,k)-mn(i-1,j-1,k))./ ...
                          (pp(uuu(ic,j,k)) + pn(uuu(i,j,k))+ ...
                          pp(vvv(i,jc,k)) + pn(vvv(i,j,k))+ ...
                          mangle(iadz,i,j,k).*...
                          (pp(www(i,j,kc)) + pn(www(i,j,k)))+eps);
            % end % for ii = 0:nx+1
        % end % for jj = 0:ny+1
    % end % for kk = 1:nzm
    if doprint == 1
        printmx(mx,ip,jp,kp);
        printmn(mn,ip,jp,kp);
    end

    % Indexing base uvw, shift (-1,-1,0) for mx,mn
    % for kk = 1:nzm
    kk = 1:nzm;
        k = kk;
        % for jj = 1:ny
        jj = 1:ny;
            j = jj+2;
            % for ii = 1:nx+1
            ii = 1:nx+1;
                i = ii+2;
                ib = i-1;
                uuu(i,j,k)=pp(uuu(i,j,k)).*...
                            min(1,min(mx(i-1,j-1,k),mn(ib-1,j-1,k))) ...
                            - pn(uuu(i,j,k)).*...
                            min(1,min(mx(ib-1,j-1,k),mn(i-1,j-1,k)));
            % end % for ii = 1:nx+1
        % end % for jj = 1:ny
    % end % for kk = 1:nzm
    if doprint == 1
        printuuu(uuu,ip,jp,kp);
    end

    % Indexing base uvw, shift (-1,-1,0) for mx,mn
    % for kk = 1:nzm
    kk = 1:nzm;
        k = kk;
        % for jj = 1:ny+1
        jj = 1:ny+1;
            j = jj+2;
            jb = j-1;
            % for ii = 1:nx
            ii = 1:nx;
                i = ii+2;
                vvv(i,j,k)=pp(vvv(i,j,k)).*...
                           min(1,min(mx(i-1,j-1,k),mn(i-1,jb-1,k))) ...
                           - pn(vvv(i,j,k)).*...
                           min(1,min(mx(i-1,jb-1,k),mn(i-1,j-1,k)));
            % end % for ii = 1:nx
        % end % for jj = 1:ny+1
    % end % for kk = 1:nzm
    if doprint == 1
        printvvv(vvv,ip,jp,kp);
    end
    
    % Indexing base uvw, shift (-1,-1,0) for mx,mn
    % for kk = 1:nzm
    kk = 1:nzm;
        k = kk;
        kb = max(1,k-1);
        % for jj = 1:ny
        jj = 1:ny;
            j = jj+2;
            % for ii = 1:nx
            ii = 1:nx;
                i = ii+2;
                www(i,j,k)=pp(www(i,j,k)).*...
                            min(1,min(mx(i-1,j-1,k), mn(i-1,j-1,kb))) ...
                             - pn(www(i,j,k)).* ...
                             min(1,min(mx(i-1,j-1,kb),mn(i-1,j-1,k)));
            % end % for ii = 1:nx
        % end % for jj = 1:ny
    % end % for kk = 1:nzm
    
    if doprint == 1
        printwww(www,ip,jp,kp);
    end

    % endif (nonos)

    %===
    % LINE 285
    %===
    % Indexing base uvw, shift (+1,+1,0) for f
    ax = zeros(nx,ny,nzm);
    ay = zeros(nx,ny,nzm);
    az = zeros(nx,ny,nzm);
    % for kk = 1:nzm
    kk = 1:nzm;
        k = kk;
        kc = k+1;
        % for jj = 1:ny
        jj = 1:ny;
            j = jj+2;
            % for  ii = 1:nx
            ii = 1:nx;
                i = ii+2;
                f(i+1,j+1,k) = max(0, f(i+1,j+1,k) - ...
                    (uuu(i+1,j,k)-uuu(i,j,k)+vvv(i,j+1,k)-vvv(i,j,k) ...
                    +(www(i,j,k+1)-www(i,j,k)).*...
                    mangle(iadz,i,j,k)).*mangle(irho,i,j,k));
            % end % for ii = 1:nx
        % end % for jj = 1:ny
    % end % for kk = 1:nzm
    
    if doprint == 1
        printf(f,ip,jp,kp);
        fprintf(1, 'rho(10) = %.6f\n', grid.rho(10));
        fprintf(1, 'rhow(10) = %.6f\n', grid.rhoi(10));
        fprintf(1, 'adz(10) = %.6f\n', grid.adz(10));
    end

    % Set ax to new field
    ax = f(4:nx+3,4:ny+3,:);
    ax = permute(ax, [3 2 1]);
    ay = NaN*ones(size(ax));
    az = NaN*ones(size(az));

end

function printmx(mx,i,j,k)
    fprintf(1, 'mx(%d,%d,%d) = %.6f\n', i,j,k, mx(i+1,j+1,k));
end
function printmn(mn,i,j,k)
    fprintf(1, 'mn(%d,%d,%d) = %.6f\n', i,j,k, mn(i+1,j+1,k));
end
function printf(f,i,j,k)
    fprintf(1, 'f(%d,%d,%d) = %.6f\n', i,j,k, f(i+3,j+3,k));
end
function printuuu(uuu,i,j,k)
    fprintf(1, 'uuu(%d,%d,%d) = %.6f\n', i,j,k, uuu(i+2,j+2,k));
end
function printvvv(vvv,i,j,k)
    fprintf(1, 'vvv(%d,%d,%d) = %.6f\n', i,j,k, vvv(i+2,j+2,k));
end
function printwww(www,i,j,k)
    fprintf(1, 'www(%d,%d,%d) = %.6f\n', i,j,k, www(i+2,j+2,k));
end