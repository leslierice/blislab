% ---------------------------------------------------------
% Plotting
% ---------------------------------------------------------
figure;

set( gcf, 'PaperSize', [3 3]);
set( gcf, 'PaperPosition', [0.25 0.25 3 3] );
set( gcf, 'Position', [0 0 420 400]);

hold;

plot( run_step6_abc_k_mt( :, 3 ), run_step6_abc_k_mt( :, 5 ), '-.', 'LineWidth', 2, 'Color',  'k');

plot( run_step6_abc_k_mt( :, 3 ), run_step6_abc_k_mt( :, 4), '.-', 'LineWidth', 1.3, 'Color',  'b' );
plot( run_step6_ab_k_mt( :, 3 ), run_step6_ab_k_mt( :, 4), '.-', 'LineWidth', 1.3, 'Color',  'r' );
plot( run_step6_naive_k_mt( :, 3 ), run_step6_naive_k_mt( :, 4), '.-', 'LineWidth', 1.3, 'Color',  'g' );


xlabel( 'k' );
ylabel( 'Effective GFLOPS (2\cdotm\cdotn\cdotk/time)' );
title( 'm=n=12000 k varies, 8 core' );

grid on;
axis square;
axis( [ 0 12000 75 200 ] );


ax = gca;
ax.YTick = [  75, 105, 135, 165, 180, 195 ];

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
print('step6_k_mt','-djpeg','-r1500')    
