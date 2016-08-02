% ---------------------------------------------------------
% Plotting
% ---------------------------------------------------------
figure;

set( gcf, 'PaperSize', [3 3]);
set( gcf, 'PaperPosition', [0.25 0.25 3 3] );
set( gcf, 'Position', [0 0 420 400]);

hold;

plot( run_step7_abc_mt( :, 1 ), run_step7_abc_mt( :, 5 ), '-.', 'LineWidth', 2, 'Color',  'k');

plot( run_step7_abc_mt( :, 1 ), run_step7_abc_mt( :, 4), '.-', 'LineWidth', 1.3, 'Color',  'b' );
plot( run_step7_ab_mt( :, 1 ), run_step7_ab_mt( :, 4), '.-', 'LineWidth', 1.3, 'Color',  'r' );
plot( run_step7_naive_mt( :, 1 ), run_step7_naive_mt( :, 4), '.-', 'LineWidth', 1.3, 'Color',  'g' );


xlabel( 'm=k=n' );
ylabel( 'Effective GFLOPS (2\cdotm\cdotn\cdotk/time)' );
title( 'm=k=n, 8 core' );

grid on;
axis square;
axis( [ 0 12000 75 195 ] );


ax = gca;
ax.YTick = [  75, 105, 135, 165, 180, 195 ];

ax.XTick = [ 0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000];
ax.XAxis.Exponent = 3;


set( gca,'FontSize', 12 );

legend( 'BLIS DSYMM', ...
        'One-level ABC Strassen', ...
        'One-level AB Strassen', ...
        'One-level Naive Strassen', ...
        'Location','SouthEast');
    
fig = gcf;
fig.PaperPositionMode = 'auto';
print('step7_mt','-djpeg','-r1500')    
