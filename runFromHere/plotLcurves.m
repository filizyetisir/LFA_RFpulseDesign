
clear all
close all


load lcurve_2spMLS_LSAR
lcurve2sp_sar = lcurve;

load lcurve_2spMLS_RF
lcurve2sp_rf = lcurve;

load lcurve_BCavg_2sp
lcurveBCavg_2sp = lcurve;

load lcurve_1spMLS_LSAR
lcurve1sp_sar = lcurve;

load lcurve_1spMLS_RF
lcurve1sp_rf = lcurve;

load lcurve_BCavg_1sp
lcurveBCavg_1sp = lcurve;

rmmse_ind = 12;
lsar_ind = 6;
gsar_ind = 7;
avgp_ind = 8;
maxp_ind = 10;
SF = 100;
DC = 0.20

figure; 

%% Local SAR
x_ind = lsar_ind; 
subplot(1,4,1);
plot(lcurveBCavg_1sp(:,rmmse_ind)*SF, lcurveBCavg_1sp(:,x_ind)*DC,'bd','MarkerFaceColor','b'); hold on; axis square; grid on; title('Local SAR')
plot(lcurveBCavg_2sp(:,rmmse_ind)*SF, lcurveBCavg_2sp(:,x_ind)*DC,'rd', 'MarkerFaceColor','r');
plot(lcurve1sp_sar(:,rmmse_ind)*SF, lcurve1sp_sar(:,x_ind)*DC,'bo-','MarkerFaceColor','b'); 
plot(lcurve1sp_rf(:,rmmse_ind)*SF, lcurve1sp_rf(:,x_ind)*DC,'bo-'); 
plot(lcurve2sp_sar(:,rmmse_ind)*SF, lcurve2sp_sar(:,x_ind)*DC,'ro-','MarkerFaceColor','r');
plot(lcurve2sp_rf(:,rmmse_ind)*SF, lcurve2sp_rf(:,x_ind)*DC,'ro-');
axis([0 40 0 67*DC])
legend('BCavg 1sp', 'BCavg 2sp', '1sp LSAR', '1sp RF', '2sp LSAR', '2sp RF')

%% Ave SAR
x_ind = gsar_ind;  
subplot(1,4,2);
plot(lcurveBCavg_1sp(:,rmmse_ind)*SF, lcurveBCavg_1sp(:,x_ind)*DC,'bd','MarkerFaceColor','b'); hold on; axis square; grid on; title('Head avg SAR')
plot(lcurveBCavg_2sp(:,rmmse_ind)*SF, lcurveBCavg_2sp(:,x_ind)*DC,'rd', 'MarkerFaceColor','r');
plot(lcurve1sp_sar(:,rmmse_ind)*SF, lcurve1sp_sar(:,x_ind)*DC,'bo-','MarkerFaceColor','b'); 
plot(lcurve1sp_rf(:,rmmse_ind)*SF, lcurve1sp_rf(:,x_ind)*DC,'bo-'); 
plot(lcurve2sp_sar(:,rmmse_ind)*SF, lcurve2sp_sar(:,x_ind)*DC,'ro-','MarkerFaceColor','r');
plot(lcurve2sp_rf(:,rmmse_ind)*SF, lcurve2sp_rf(:,x_ind)*DC,'ro-');
axis([0 40 0 22*DC])


%% AVG power
x_ind = avgp_ind; 
subplot(1,4,3); 
plot(lcurveBCavg_1sp(:,rmmse_ind)*SF, lcurveBCavg_1sp(:,x_ind)*DC,'bd','MarkerFaceColor','b'); hold on; axis square; grid on; title('Max avg power per channel')
plot(lcurveBCavg_2sp(:,rmmse_ind)*SF, lcurveBCavg_2sp(:,x_ind)*DC,'rd', 'MarkerFaceColor','r');
plot(lcurve1sp_sar(:,rmmse_ind)*SF, lcurve1sp_sar(:,x_ind)*DC,'bo-','MarkerFaceColor','b'); 
plot(lcurve1sp_rf(:,rmmse_ind)*SF, lcurve1sp_rf(:,x_ind)*DC,'bo-'); 
plot(lcurve2sp_sar(:,rmmse_ind)*SF, lcurve2sp_sar(:,x_ind)*DC,'ro-','MarkerFaceColor','r');
plot(lcurve2sp_rf(:,rmmse_ind)*SF, lcurve2sp_rf(:,x_ind)*DC,'ro-');
axis([0 40 0 67*DC])


%% MAX power
x_ind = maxp_ind;
subplot(1,4,4); 
plot(lcurveBCavg_1sp(:,rmmse_ind)*SF, lcurveBCavg_1sp(:,x_ind),'bd','MarkerFaceColor','b'); hold on; axis square; grid on; title('Max peak power per channel')
plot(lcurveBCavg_2sp(:,rmmse_ind)*SF, lcurveBCavg_2sp(:,x_ind),'rd', 'MarkerFaceColor','r');
plot(lcurve1sp_sar(:,rmmse_ind)*SF, lcurve1sp_sar(:,x_ind),'bo-','MarkerFaceColor','b'); 
plot(lcurve1sp_rf(:,rmmse_ind)*SF, lcurve1sp_rf(:,x_ind),'bo-'); 
plot(lcurve2sp_sar(:,rmmse_ind)*SF, lcurve2sp_sar(:,x_ind),'ro-','MarkerFaceColor','r');
plot(lcurve2sp_rf(:,rmmse_ind)*SF, lcurve2sp_rf(:,x_ind),'ro-');
axis([0 40 0 225])

