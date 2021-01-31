%--------------------------------------------------------------------------
% Script to generate Figure 2
% 1) Load repeated evaluations of Log Likelihoods
% 2) Plot Results
%--------------------------------------------------------------------------

%% Housekeeping
clear; clc;
close all
set(0,'defaultTextInterpreter','latex');

%% Load Results

nFilter1  = 'COPF400';
nFilter2  = 'BSPF1000';
ME_scale  = '15';  % 10 or 15
nDataSet  = 'Sim3';

lhFolder = sprintf(strcat('ME',ME_scale,'_',nDataSet));
lhPath   = [pwd, '\LikelihoodEvaluations\', lhFolder ,'\'];

% row1: standard deviation of loglh
% row2: average of loglh
mLogStat_Filter1 = csvread( strcat(lhPath, nFilter1,'_LogStats.csv') );
mLogStat_Filter2 = csvread( strcat(lhPath, nFilter2,'_LogStats.csv') );

vMeanFilter1  = mLogStat_Filter1(2,:);
vSDFilter1    = mLogStat_Filter1(1,:);
vSDFilter2    = mLogStat_Filter2(1,:);

%% Generate the plot

figure(1);
clf;
set(figure(1),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);

hold on
scatter(vMeanFilter1,vSDFilter1,100,'b','filled')
scatter(vMeanFilter1,vSDFilter2,150,'r','*' )
xlabel('Mean(ln p(Y|\theta^i))', 'Interpreter','tex')
ylabel('StdD(ln p(Y|\theta^i))', 'Interpreter','tex')
set(gca,'FontSize',30)
box on

% For ME10
%axis([-340 -290 0 10]) 

% For ME15
axis([-360 -320 0 10]) 
legend('COPF M=400','BSPF M=1,000','Location','northwest')

saveas(figure(1), [lhPath '\LhSD_Scatter.pdf'] );

