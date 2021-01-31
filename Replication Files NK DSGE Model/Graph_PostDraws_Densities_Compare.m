%% Script to compare posterior densities

clc
clear all;
close all;


nFilter1      = 'COPF'; % Choose BSPF, COPF, COPFexactZLB or KF
nMhrun1       = '1';
nPrior1       = '1';

nFilter2      = 'BSPF';
nMhrun2       = '1';
nPrior2       = '1';

nData   = 'Sim3'; % Sim1 or US


% Change the path for desired draws
postPath1 = [pwd, strcat('\PosteriorDraws\',nFilter1,'_Prior',nPrior1,'_',nData,'_Mhrun',nMhrun1,'\')];
postPath2 = [pwd, strcat('\PosteriorDraws\',nFilter2,'_Prior',nPrior2,'_',nData,'_Mhrun',nMhrun2,'\')];
figPath   = [pwd, strcat('\PosteriorDraws\Figures_',nFilter1,'Mhrun',nMhrun1,'_',nFilter2,'Mhrun',nMhrun2,'_Prior',nPrior1,'_',nData,'\')];
[~, ~, ~] = mkdir(figPath);

% Read parameter draws and likelihood values
mParam_Candidate1 = csvread([postPath1, strcat(nFilter1,'_Prior',nPrior1,'_',nData,'_Mhrun',nMhrun1,'_CandidateDraws.csv')]);
mParam_Accepted1  = csvread([postPath1, strcat(nFilter1,'_Prior',nPrior1,'_',nData,'_Mhrun',nMhrun1,'_Draws.csv')]);
vLogLH_Accepted1  = csvread([postPath1, strcat(nFilter1,'_Prior',nPrior1,'_',nData,'_Mhrun',nMhrun1,'_LogPosterior.csv')]);
vLogLH_Accepted1  = vLogLH_Accepted1(:,1);

mParam_Accepted2  = csvread([postPath2, strcat(nFilter2,'_Prior',nPrior2,'_',nData,'_Mhrun',nMhrun2,'_Draws.csv')]);
vLogLH_Accepted2  = csvread([postPath2, strcat(nFilter2,'_Prior',nPrior2,'_',nData,'_Mhrun',nMhrun2,'_LogPosterior.csv')]);
vLogLH_Accepted2  = vLogLH_Accepted2(:,1);

% Read true parameters
%vTrueParam = csvread([pwd, '\PosteriorDraws\BSPF_Prior1_Sim1_Mhrun1_InitialValue.csv']);
vTrueParam = csvread([pwd, '\Data\COPFexactZLB_DGP1_Param.csv']);
vTrueParam = vTrueParam';
                                                                                    
%% Plot densities of parameters

sNameParam = ["tau", "kappa", "psi1", "psi2", "rho_R", "rho_g", "rho_d", "rho_z", "sig_r", ...
              "sig_g", "sig_d", "sig_z", "eta", "nu", "chi_h", "gstar", "rAnet", "gamQnet", ...
              "piAnet", "pibarAnet", "er0", "g0", "z0", "d0", "y0", "pi0", "c0", "R0", "ylag0"];
nParamID   = size(sNameParam,2);
vTrueParam = mParam_Candidate1(1,:); % only used for Sim1 data set, sampler is initialized at true para

vFix = zeros(nParamID,1);
for iParam = 1:nParamID
    if mParam_Candidate1(2,iParam) == mParam_Candidate1(3,iParam)
        vFix(iParam,1) = 1;
    end
end

for iParam = 1:nParamID
    
    if vFix(iParam,1) == 0
        
        set(figure(1),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
        [vY,vX,bw] = ksdensity(mParam_Accepted1(:,iParam));
        bw = 3*bw;
        [vY,vX] = ksdensity(mParam_Accepted1(:,iParam),'Bandwidth',bw);
        plot(vX,vY,'Color','b','LineWidth',4)
        vY1 = vY;
        hold on
        [vY,vX,bw] = ksdensity(mParam_Accepted2(:,iParam));
        bw = 4*bw;
        [vY,vX] = ksdensity(mParam_Accepted2(:,iParam),'Bandwidth',bw);        
        plot(vX,vY,'--','Color','r','LineWidth',4)
        hold on
        if strcmp(nData,'US') == 0
            line([vTrueParam(iParam) vTrueParam(iParam)], [0 max([vY1,vY])*1.1],'Color','k','LineWidth',2.5);
        end
        axis tight
        sTitle = sprintf(sNameParam(iParam));
        %title(sTitle)
        legend(nFilter1,nFilter2,'location','best')
        set(gca,'FontSize',30)    
        sName = sprintf('PostDensity_%s.pdf',sNameParam(iParam));
        saveas(figure(1), [figPath, sName] ); 
        close all
    
    end    
    
end

%% Plot densities of log likelihood values

set(figure(1),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
[vY,vX,bw] = ksdensity(vLogLH_Accepted1);
bw = 3*bw;
[vY,vX,bw] = ksdensity(vLogLH_Accepted1,'Bandwidth',bw); 
plot(vX,vY,'Color','b','LineWidth',2)
hold on
[vY,vX,bw] = ksdensity(vLogLH_Accepted2);
bw = 4*bw;
[vY,vX,bw] = ksdensity(vLogLH_Accepted2,'Bandwidth',bw); 
plot(vX,vY,'--','Color','r','LineWidth',2)
axis tight
title("Density of Log Likelihood Values")
legend(nFilter1,nFilter2)
set(gca,'FontSize',20)    
saveas(figure(1), [figPath, 'AccLogLH.pdf'] ); 
close all

