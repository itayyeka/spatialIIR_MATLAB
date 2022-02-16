function beautifyFigureAndMaximize
axes = gca;
fig = gcf;
h = findobj ('type', 'line');
set (h(1:end), 'LineWidth', 2)
set (get (axes, 'XLabel'), 'FontSize', 20);
set (get (axes, 'YLabel'), 'FontSize', 20);
set (get (axes, 'Title'), 'FontSize', 24);
set (axes, 'Box', 'off', 'TickDir', 'out', 'XMinorTick', 'on', 'YMinorTick', 'off', 'FontSize', 16)
screen_size = get(0, 'ScreenSize');
set(fig, 'Position', [0 0 screen_size(3) screen_size(4)]);
end