% ---------------------------------------------------------
% Plotting
% ---------------------------------------------------------
figure;

set( gcf, 'PaperSize', [550 550]);
set( gcf, 'PaperPosition', [0.25 0.25 3 3] );
set( gcf, 'Position', [0 0 420 400]);

hold;

plot( run_step6_abc_st( :, 1 ), run_step6_abc_st( :, 5 ), '-.', 'LineWidth', 2, 'Color',  'k');

plot( run_step6_abc_st( :, 1 ), run_step6_abc_st( :, 4), '.-', 'LineWidth', 1.3, 'Color',  'b' );
plot( run_step6_ab_st( :, 1 ), run_step6_ab_st( :, 4), '.-', 'LineWidth', 1.3, 'Color',  'r' );
plot( run_step6_naive_st( :, 1 ), run_step6_naive_st( :, 4), '.-', 'LineWidth', 1.3, 'Color',  'g' );


xlabel( 'm=k=n' );
ylabel( 'Effective GFLOPS (2\cdotm\cdotn\cdotk/time)' );
title( 'm=k=n, 1 core' );

grid on;
axis square;
axis( [ 0 12000 18 26 ] );

ax = gca;
ax.YTick = [  18, 19 20, 21, 22, 23, 23.76, 24, 25, 26 ];

ax.XTick = [ 0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000];
ax.XAxis.Exponent = 3;


set( gca,'FontSize', 12 );

legend( 'BLIS DGEMM', ...
        'One-level ABC Strassen', ...
        'One-level AB Strassen', ...
        'One-level Naive Strassen', ...
        'Location','SouthEast');
    
fig = gcf;
fig.PaperPositionMode = 'auto';
print('step6_st','-djpeg','-r1500')
