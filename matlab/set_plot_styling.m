%% Set plot styling
set(0,'DefaultLegendAutoUpdate','off')
if latex ~= 0
    set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
    set(0, 'DefaultLegendInterpreter', 'latex');
    set(0, 'DefaultTextInterpreter', 'latex');
else
    set(0, 'DefaultAxesTickLabelInterpreter', 'default');
    set(0, 'DefaultLegendInterpreter', 'default');
    set(0, 'DefaultTextInterpreter', 'default');
end