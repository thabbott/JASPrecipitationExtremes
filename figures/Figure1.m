% ABC2020 Figure 1

% Setup
addpath(genpath('..\matlab'));
SAM_defineConstants;
SOM_defineConstants;

%% Create profiles
w1 = @(p,ps) -sin(pi*p/ps); % First baroclinic mode
w2 = @(p,ps) -sin(2*pi*p/ps); % Second baroclinic modes
ps = 100000; % Surface pressure, Pa
p = linspace(ps, 1, 10000); % Pressure profile, Pa
wp1 = w1(p, ps); % omega profile 1
wp2 = (w1(p,ps) + 0.5*w2(p,ps)); wp2 = -wp2/min(wp2); % omega profile 2
wp3 = (w1(p,ps) - 0.5*w2(p,ps)); wp3 = -wp3/min(wp3); % omega profile 3

%% Compute scaling (w/ liquid CC thermodynamics)
SSTs = linspace(280, 310, 200);
qv0 = SAM_qsatWater(SSTs, ps);
wrangle = @(x) repmat(reshape(x, [1 1 numel(x)]), [size(SSTs) 1]);
C1 = CC_condensateIntegral(SSTs, qv0, wrangle(wp1),...
    wrangle(p), 1);
C2 = CC_condensateIntegral(SSTs, qv0, wrangle(wp2),...
    wrangle(p), 1);
C3 = CC_condensateIntegral(SSTs, qv0, wrangle(wp3),...
    wrangle(p), 1);

pp = reshape(p, [1 1 numel(p)]);
[T1,~,dpqv1] = CC_pseudoAdiabat(280, SAM_qsatWater(280,ps), pp, 1);
[T2,~,dpqv2] = CC_pseudoAdiabat(290, SAM_qsatWater(290,ps), pp, 1);
[T3,~,dpqv3] = CC_pseudoAdiabat(300, SAM_qsatWater(300,ps), pp, 1);
[T4,~,dpqv4] = CC_pseudoAdiabat(310, SAM_qsatWater(310,ps), pp, 1);
T1 = squeeze(T1); T2 = squeeze(T2); T3 = squeeze(T3); T4 = squeeze(T4);
dpqv1 = squeeze(dpqv1); dpqv2 = squeeze(dpqv2); dpqv3 = squeeze(dpqv3);
dpqv4 = squeeze(dpqv4);

%% Plot results

%% Mass flux
figure(standard_figure('quarter-page')); hold on;
plot(-wp1, p/100, 'k', 'LineWidth', 2);
xlabel('\omega (norm.)');
ylabel('Pressure (hPa)');
set(gca, 'YDir', 'reverse');
plot(-wp2, p/100, 'b', 'LineWidth', 2);
plot(-wp3, p/100, 'r', 'LineWidth', 2);
title('(b)');
pngandpdf('Figure1B');

%% Moisture lapse rate
cols = viridis(4);
figure(standard_figure('quarter-page')); hold on;
plot(1e5*dpqv1, p/100, 'LineWidth', 2, 'Color', cols(1,:));
xlabel('dq*/dp (g/kg/hPa)');
ylabel('Pressure (hPa)');
set(gca, 'YDir', 'reverse');
plot(1e5*dpqv2, p/100, 'b', 'LineWidth', 2, 'Color', cols(2,:));
plot(1e5*dpqv3, p/100, 'r', 'LineWidth', 2, 'Color', cols(3,:));
plot(1e5*dpqv4, p/100, 'r', 'LineWidth', 2, 'Color', cols(4,:));
text(0.06, 750, '280 K', 'Color', cols(1,:), 'FontSize', 10);
text(0.06, 815, '290 K', 'Color', cols(2,:), 'FontSize', 10);
text(0.06, 880, '300 K', 'Color', cols(3,:), 'FontSize', 10);
text(0.06, 945, '310 K', 'Color', cols(4,:), 'FontSize', 10);
title('(a)');
pngandpdf('Figure1A');

%% Scaling
figure(standard_figure('half-page')); hold on;
scal = @(x,y) 100*2*(x - y)./(x + y)./diff(SSTs);
sstm = 0.5*(SSTs(2:end) + SSTs(1:end-1));
plot(sstm, 100*SAM_L_c./(SAM_R_v.*sstm.^2), 'k--', 'LineWidth', 1.5);
ylim([0 12]);
xticks([280 285 290 295 300 305 310]);
xlabel('SST (K)');
ylabel('Thermodynamic mode (% K^{-1})');
l = [];
l = [l, ...
    plot(sstm, scal(C1(2:end), C1(1:end-1)), 'k', 'LineWidth', 2)];
N = numel(C2);
l = [l, ...
    plot(sstm, scal(C2(2:N), C2(1:N-1)), 'b', 'LineWidth', 2)];
l = [l, ...
    plot(sstm, scal(C3(2:N), C3(1:N-1)), 'r', 'LineWidth', 2)];
title('(c)');
leg = legend(l([2 1 3]), {'Top-heavy', 'Neutral', 'Bottom-heavy'}, ...
    'Location', 'southwest');
leg.Box = 'off';
leg.ItemTokenSize = [10, 1];
pngandpdf('Figure1C');
