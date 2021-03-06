%% Generate figures with standard screen-independent sizes
function fig = standard_figure(ftype)

    fig = figure('Units', 'Inches', 'PaperUnits', 'Inches');
    switch ftype
        case 'medium-very-narrow'
            set(fig, 'Position', [1 1 4 6], 'PaperSize', [8 6]);
        case 'medium'
            set(fig, 'Position', [1 1 8 6], 'PaperSize', [8 6]);
        case 'medium-wide'
            set(fig, 'Position', [1 1 12 6], 'PaperSize', [12 6]);
        case 'medium-narrow'
            set(fig, 'Position', [1 1 6 6], 'PaperSize', [6 6]);
        case 'wide'
            set(fig, 'Position', [1 1 12 6], 'PaperSize', [12 6]);
        case 'half-page'
            set(fig, 'Position', [1 1 3.75 2.82], 'PaperSize', [3.75 2.82]);
        case 'half-page-short'
            set(fig, 'Position', [1 1 3.75 1.75], 'PaperSize', [3.75 1.75]);
        case 'quarter-page'
            set(fig, 'Position', [1 1 1.875 2.82], 'PaperSize', [1.875 2.82]);
        case 'third-page-tall'
            set(fig, 'Position', [1 1 2.25 5.64], 'PaperSize', [2.25 5.64]);
        case 'third-page'
            set(fig, 'Position', [1 1 2.25 2.82], 'PaperSize', [2.25 2.82]);
        otherwise
            error('Unrecognized standard figure: %s', ftype);
    end

end
