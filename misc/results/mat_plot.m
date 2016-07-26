%clear;

%step6;

% ---------------------------------------------------------
% Plotting
% ---------------------------------------------------------
figure;
% hFig = figure(1);

set( gcf, 'PaperSize', [3 3]);
set( gcf, 'PaperPosition', [0.25 0.25 3 3] );
set( gcf, 'Position', [0 0 420 400]);

hold;

plot( run_step6_st( :, 1 ), run_step6_st( :, 5 ), '-.', 'LineWidth', 2, 'Color',  'k');
%plot( run_step7_st( :, 1 ), run_step7_st( :, 5 ), '-.', 'LineWidth', 2, 'Color',  'k');
%plot( run_step6_mt( :, 1 ), run_step6_mt( :, 5 ), '-.', 'LineWidth', 2, 'Color',  'k');
%plot( run_step7_mt( :, 1 ), run_step7_mt( :, 5 ), '-.', 'LineWidth', 2, 'Color',  'k');
%plot( sb_mkl_gemm( :, 1 ), sb_mkl_gemm( :, 5 ), '.-', 'LineWidth', 1.3, 'Color',  'k');

plot( run_step6_st( :, 1 ), run_step6_st( :, 4), '.-', 'LineWidth', 1.3, 'Color',  'b' );
%plot( run_step7_st( :, 1 ), run_step7_st( :, 4), '.-', 'LineWidth', 1.3, 'Color',  'b' );
%plot( run_step6_mt( :, 1 ), run_step6_mt( :, 4), '.-', 'LineWidth', 1.3, 'Color',  'b' );
%plot( run_step7_mt( :, 1 ), run_step7_mt( :, 4), '.-', 'LineWidth', 1.3, 'Color',  'b' );

%plot( sb_stra_1level2( :, 1 ), sb_stra_1level2( :, 8), '.-', 'LineWidth', 3, 'Color',  'b' );
%plot( sb_stra_1levelref_par( :, 1 ), sb_stra_1levelref_par( :, 8), '--', 'LineWidth', 1.3, 'Color',  [0 0.2 1.0] );


xlabel( 'm=k=n' );
ylabel( 'Effective GFLOPS (2\cdotm\cdotn\cdotk/time)' );
%title( 'm=k=n, 1 core' );
title( 'm=k=n, 8 core' );

grid on;
axis square;
axis( [ 0 12000 18 26 ] );
%axis( [ 0 12000 125 195 ] );


ax = gca;
ax.YTick = [  18, 19 20, 21, 22, 23, 23.76, 24, 25, 26 ];
%ax.YTick = [  125, 165, 180, 195 ];

ax.XTick = [ 0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000];
ax.XAxis.Exponent = 3;


set( gca,'FontSize', 12 );

legend( 'BLIS DGEMM', ...
        'One-level ABC Strassen', ...
        'Location','SouthEast');


