%% ABC2020 Figure A1

% Setup
addpath(genpath('..\matlab'));
SAM_defineConstants;
SOM_defineConstants;

%% Create profiles
ps = 100000; % Surface pressure, Pa
p = linspace(ps, 100, 10000); % Pressure profile, Pa
SSTs = linspace(280, 310, 200);
qv0 = SAM_qsatWater(SSTs, ps);
wrangle = @(x) repmat(reshape(x, [1 1 numel(x)]), [size(SSTs) 1]);
pp = reshape(p, [1 1 numel(p)]);
[T1,~,dpqv1] = CC_pseudoAdiabat(280, SAM_qsatWater(280,ps), pp, 1);
[T3,~,dpqv3] = CC_pseudoAdiabat(295, SAM_qsatWater(295,ps), pp, 1);
[T4,~,dpqv4] = CC_pseudoAdiabat(310, SAM_qsatWater(310,ps), pp, 1);
T1 = squeeze(T1); T3 = squeeze(T3); T4 = squeeze(T4);
dpqv1 = squeeze(dpqv1); dpqv3 = squeeze(dpqv3);
dpqv4 = squeeze(dpqv4);
e1 = SAM_psatWater(T1);
i = find(e1 == min(e1), 1);
e1(i+1:end) = 0;
q1 = e1./p';
qp1 = q1./p';
phiT1 = ((SAM_L_c * SAM_R / (SAM_R_v * SAM_c_p)) ./ T1) .* ...
    ((1 + (SAM_L_c.*q1)./(SAM_R.*T1)) ./ ...
    (1 + q1.*(SAM_c_pv./SAM_c_p + (SAM_L_c./(SAM_c_p.*T1)).*...
    ((SAM_L_c./(SAM_R_v.*T1)) - 1)))) - 1;
e3 = SAM_psatWater(T3);
i = find(e3 == min(e3), 1);
e3(i+1:end) = 0;
q3 = e3./p';
qp3 = q3./p';
phiT3 = ((SAM_L_c * SAM_R / (SAM_R_v * SAM_c_p)) ./ T3) .* ...
    ((1 + (SAM_L_c.*q3)./(SAM_R.*T3)) ./ ...
    (1 + q3.*(SAM_c_pv./SAM_c_p + (SAM_L_c./(SAM_c_p.*T3)).*...
    ((SAM_L_c./(SAM_R_v.*T3)) - 1)))) - 1;
e4 = SAM_psatWater(T4);
i = find(e4 == min(e4), 1);
e4(i+1:end) = 0;
q4 = e4./p';
qp4 = q4./p';
phiT4 = ((SAM_L_c * SAM_R / (SAM_R_v * SAM_c_p)) ./ T4) .* ...
    ((1 + (SAM_L_c.*q4)./(SAM_R.*T4)) ./ ...
    (1 + q4.*(SAM_c_pv./SAM_c_p + (SAM_L_c./(SAM_c_p.*T4)).*...
    ((SAM_L_c./(SAM_R_v.*T4)) - 1)))) - 1;

%% Plot results
cols = [0 0 0; 0 0 1];
figure(standard_figure('half-page')); hold on;
l = [plot(1e5*dpqv1, p/100, 'LineWidth', 2, 'Color', cols(1,:))];
plot(1e5*SAM_R/SAM_R_v*qp1, p/100, '-', ...
    'LineWidth', 1, 'Color', cols(1,:));
plot(1e5*SAM_R/SAM_R_v*e1/p(1)^2, p/100, '--', 'LineWidth', 1, 'Color', cols(1,:));
xlabel('dq*/dp (g/kg/hPa)');
ylabel('Pressure (hPa)');
set(gca, 'YDir', 'reverse');
xlim([-0.0005 0.08]);
% l = [l, plot(1e5*dpqv3, p/100, 'r', 'LineWidth', 2, 'Color', cols(2,:))];
% plot(1e5*SAM_R/SAM_R_v*qp3, p/100, '-', 'LineWidth', 1, 'Color', cols(2,:));
% plot(1e5*SAM_R/SAM_R_v*e3/p(1)^2, p/100, '--', 'LineWidth', 1, 'Color', cols(2,:));
l = [l, plot(1e5*dpqv4, p/100, 'r', 'LineWidth', 2, 'Color', cols(2,:))];
plot(1e5*SAM_R/SAM_R_v*qp4, p/100, '-', 'LineWidth', 1, 'Color', cols(2,:));
plot(1e5*SAM_R/SAM_R_v*e4/p(1)^2, p/100, '--', 'LineWidth', 1, 'Color', cols(2,:));
l = [l, plot([nan, nan], [nan, nan], 'k-')];
l = [l, plot([nan, nan], [nan, nan], 'k--')];
leg = legend(l, {'280 K', '310 K', 'M_v e^*/(M_a p^2)', 'e^* (rescaled)'}, ...
    'Position', [0.5972    0.1746    0.3130    0.3007]);
leg.Box = 'off';
leg.ItemTokenSize = [20, 1];
title('(a)');
pngandpdf('FigureA1a');
%% 
figure(standard_figure('half-page')); hold on;
plot(phiT1, p/100, 'LineWidth', 2, 'Color', cols(1,:));
set(gca, 'YDir', 'reverse');
xlim([0, 10]);
% plot(phiT3, p/100, 'LineWidth', 2, 'Color', cols(2,:));
plot(phiT4, p/100, 'LineWidth', 2, 'Color', cols(2,:));
ylabel('Pressure (hPa)');
xlabel('(T_0/T)\phi(T,p) - 1 (nondim.)');
leg = legend({'280 K', '310 K'}, 'Location', 'southeast');
leg.Box = 'off';
title('(b)');
leg.ItemTokenSize = [20, 1];
pngandpdf('FigureA1b');
