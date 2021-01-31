%% Script to plot ACFs of Posterior Draws

clc
clear all;
close all;

%% Read posterior draws 

nFilter1      = 'COPF'; % Choose BSPF, COPF, COPFexactZLB or KF
nMhrun1       = '1';

nFilter2      = 'BSPF';
nMhrun2       = '1';

nData   = 'Sim3'; % Sim1 or US
nPrior  = '1';

% Set the path for posterior draws and figures
postPath1 = [pwd, strcat('\PosteriorDraws\',nFilter1,'_Prior',nPrior,'_',nData,'_Mhrun',nMhrun1,'\')];
postPath2 = [pwd, strcat('\PosteriorDraws\',nFilter2,'_Prior',nPrior,'_',nData,'_Mhrun',nMhrun2,'\')];
figPath   = [pwd, strcat('\PosteriorDraws\Figures_',nFilter1,'Mhrun',nMhrun1,'_',nFilter2,'Mhrun',nMhrun2,'_Prior',nPrior,'_',nData,'\')];
[~, ~, ~] = mkdir(figPath);

mParam_Candidate1 = csvread([postPath1, strcat(nFilter1,'_Prior',nPrior,'_',nData,'_Mhrun',nMhrun1,'_CandidateDraws.csv')]);
mParam_Accepted1  = csvread([postPath1, strcat(nFilter1,'_Prior',nPrior,'_',nData,'_Mhrun',nMhrun1,'_Draws.csv')]);
vLogLH_Accepted1  = csvread([postPath1, strcat(nFilter1,'_Prior',nPrior,'_',nData,'_Mhrun',nMhrun1,'_LogPosterior.csv')]);

mParam_Accepted2  = csvread([postPath2, strcat(nFilter2,'_Prior',nPrior,'_',nData,'_Mhrun',nMhrun2,'_Draws.csv')]);
vLogLH_Accepted2  = csvread([postPath2, strcat(nFilter2,'_Prior',nPrior,'_',nData,'_Mhrun',nMhrun2,'_LogPosterior.csv')]);


%% Prepare output files and compute ACFs

sNameParam = ["tau", "kappa", "psi1", "psi2", "rho_R", "rho_g", "rho_d", ...
              "rho_z", "sig_r", "sig_g", "sig_d", "sig_z", "eta", "nu", "chi_h", ...
              "gstar", "rAnet", "gamQnet", "piAnet", "pibarAnet", "er0", ...
              "g0", "z0", "d0", "y0", "pi0", "c0", "R0", "ylag0"];
          
nParamID = size(sNameParam,2);
vTrueParam = mParam_Candidate1(1,:);

vFix = zeros(nParamID,1);
for iParam = 1:nParamID
    if mParam_Candidate1(2,iParam) == mParam_Candidate1(3,iParam)
        vFix(iParam,1) = 1;
    end
end

nMaxLag = 100;
mACF_Filter1 = zeros(nMaxLag+1,nParamID);
mACF_Filter2 = zeros(nMaxLag+1,nParamID);

for iParam = 1:nParamID 
    for iLag = 1:nMaxLag+1    
        mACF_Filter1(iLag,iParam) = corr( mParam_Accepted1(1:end-nMaxLag,iParam), mParam_Accepted1(iLag:end-nMaxLag+iLag-1,iParam) );
        mACF_Filter2(iLag,iParam) = corr( mParam_Accepted2(1:end-nMaxLag,iParam), mParam_Accepted2(iLag:end-nMaxLag+iLag-1,iParam) );
    end
    sNameParam(iParam)
    mACF_Filter1(nMaxLag+1,iParam)
    mACF_Filter2(nMaxLag+1,iParam)
end

%% Plot  ACFs of posterior draws

for iParam = 1:nParamID
    
    if vFix(iParam,1) == 0
        
        set(figure(1),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
        hold on
        plot(0:nMaxLag,[mACF_Filter1(:,iParam)],'Color','b','LineWidth',4)
        plot(0:nMaxLag,[mACF_Filter2(:,iParam)],'--','Color','r','LineWidth',4)
        legend(nFilter1,nFilter2,'Location','southwest')
        hold off
        box on
        axis([0 100 0 1])
        xlabel('Lag')
        set(gca,'FontSize',30)
%         sTitle = sprintf(sNameParam(iParam));
%         title(sTitle)
%         set(gca,'FontSize',30)    
        sName = sprintf('ACF_%s.pdf',sNameParam(iParam));
        saveas(figure(1), [figPath, sName] ); 
        close all
       
    
    end    
    
end


