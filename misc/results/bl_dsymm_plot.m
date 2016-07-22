PEAK_PERF=25.15;

% ---------------------------------------------------------
% Plotting
% ---------------------------------------------------------
figure;
% hFig = figure(1);
%set(hFig, 'Position', [0 0 160 240])

set( gcf, 'PaperSize', [33 3])
set( gcf, 'PaperPosition', [0.25 0.25 3 3] );
set( gcf, 'Position', [0 0 600 400]);

hold;

plot( run_step7_st( :, 1 ), run_step7_st( :,4), '.-', 'LineWidth', 2, 'Color',  [0 0.2 1.0] );
plot( run_step7_st( :, 1 ), run_step7_st( :, 5), '.-', 'LineWidth', 2, 'Color', [1 0 0.2] );

xlabel( 'm=k=n' );
ylabel( 'Effective GFLOPS (2*m*k/time)' );
title( 'DSYMM (m=k=n, 1 core)' );

grid on;
axis square;
axis( [ 0 12000 20 PEAK_PERF ] );

ax = gca;
ax.YTick = [  20, 21, 22, 23, 24, 24.99 ];

ax.XTick = [ 0, 2000, 4000, 6000, 8000, 10000, 12000];
set( gca,'FontSize',14 );

legend( 'One-level Strassen', ...
        'Blislab DSYMM', ...
        'Location','SouthEast');

