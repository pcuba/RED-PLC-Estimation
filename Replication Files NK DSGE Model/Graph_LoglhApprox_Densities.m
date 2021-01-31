%--------------------------------------------------------------------------
% Script to generate Figure 1
% Run separately for each choice of ME
%--------------------------------------------------------------------------

clear; 
clc;
close all

set(0,'defaultTextInterpreter','latex');

LHFolder = 'LikelihoodEvaluations\ME10_Sim3';

LHPath = [pwd, '\', LHFolder ,'\'];

dataBSPF1 = csvread( [LHPath, 'BSPF1000_LogLH.csv'] );
dataBSPF2 = csvread( [LHPath, 'BSPF10000_LogLH.csv'] );
dataCOPF  = csvread( [LHPath, 'COPF400_LogLH.csv'] );

dataBSPF1(dataBSPF1==0) = nan;
vRangeLogLH_BSPF1 = prctile(dataBSPF1, 95) -prctile(dataBSPF1, 5);
dataBSPF2(dataBSPF2==0) = nan;
vRangeLogLH_BSPF2 = prctile(dataBSPF2, 95) -prctile(dataBSPF2, 5);
dataCOPF(dataCOPF==0) = nan;
vRangeLogLH_COPF = prctile(dataCOPF, 95) -prctile(dataCOPF, 5);

% Select parameters for which to plot densities
% we only use the first parameter, which is the true value used to generate
% the data 

vPara = [ 1 ];

for iPara = 1:length(vPara)
    
    nPara = vPara(iPara);

    vLogLH_BSPF1 = dataBSPF1(:,nPara);
    vLogLH_BSPF2 = dataBSPF2(:,nPara);    
    vLogLH_COPF  = dataCOPF(:,nPara);

    [likDistY_BSPF1,likDistX_BSPF1] = ksdensity(vLogLH_BSPF1);
    [likDistY_BSPF2,likDistX_BSPF2] = ksdensity(vLogLH_BSPF2);    
    [likDistY_COPF,likDistX_COPF] = ksdensity(vLogLH_COPF);     

    figure(1);clf;
    set(figure(1),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
    plot(likDistX_BSPF1,likDistY_BSPF1,'Color','r','LineWidth',4)
    hold on
    plot(likDistX_BSPF2,likDistY_BSPF2,'LineStyle','--','Color','r','LineWidth',4)
    hold on    
    plot(likDistX_COPF,likDistY_COPF,'Color','b','LineWidth',4)

    % for ME15 & sample 401-540
    %axis([-345 -325 0 0.8]) 
    %legend('BSPF M=1,000','BSPF M=10,000','COPF M=400','Location','northwest')
    
    % for ME10 & sample 401-540
    axis([-320 -300 0 0.8]) 
        
    
    set(gca,'FontSize',30)
    %legend('Bootstrap PF','Cond.Opt. PF','Location','northwest')
    sName = sprintf('LHDist_Para%d.pdf',nPara);
    saveas(figure(1), [pwd '\' LHFolder '\' sName] );
    %close all

end






