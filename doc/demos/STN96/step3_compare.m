clearvars, clearvars -global, clc

niiNMSE_L = load_untouch_nii( fullfile('scan1','Tracking','PROB','Results_LIFE','fit_NRMSE.nii') );
niiNMSE_C = load_untouch_nii( fullfile('scan1','Tracking','PROB','Results_COMMIT','fit_NRMSE.nii') );
niiMASK   = load_untouch_nii( fullfile('scan1','Tracking','PROB','dictionary_mask.nii') );

figure(1), set(gcf,'Color',[1 1 1]);
cla, set(gca,'Color',[.97 .97 .97]); hold on
x = linspace(0,1,100);
yL = hist( niiNMSE_L.img(niiMASK.img>0), x );
yC = hist( niiNMSE_C.img(niiMASK.img>0), x );
plot( x, yL, '- ', 'LineWidth', 3, 'Color',[.8 0 0] )
plot( x, yC, '- ', 'LineWidth', 3, 'Color',[0 .8 0] )
grid on, box on
legend( 'LiFE', 'COMMIT' )
xlabel( 'NRMSE' ), ylabel('frequency')
title('LiFE vs COMMIT')

figure(2), set(gcf,'Color',[1 1 1]);
cla, set(gca,'Color',[.97 .97 .97]); hold on
yL = niiNMSE_L.img( niiMASK.img>0 );
yC = niiNMSE_C.img( niiMASK.img>0 );
plot( yL, yC, 'bx' )
plot( [0 1], [0 1], 'k--', 'LineWidth', 2 )
grid on, box on
axis([0 1 0 1])
xlabel( 'LiFE' ), ylabel( 'COMMIT' )
