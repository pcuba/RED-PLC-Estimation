%% Script to Generate Scatter Plot of ACFs

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
vNotFixed  = ones(nParamID,1);

for iParam = 1:nParamID
    if mParam_Candidate1(2,iParam) == mParam_Candidate1(3,iParam)
        vNotFixed(iParam,1) = 0;
    end
end

nMaxLag = 100;
% mACF_Filter1 = zeros(nMaxLag+1,nParamID);
% mACF_Filter2 = zeros(nMaxLag+1,nParamID);
% 
% for iParam = 1:nParamID 
%     if vNotFixed(iParam,1) == 1
%        for iLag = 1:nMaxLag+1    
%            mACF_Filter1(iLag,iParam) = corr( mParam_Accepted1(1:end-nMaxLag,iParam), mParam_Accepted1(iLag:end-nMaxLag+iLag-1,iParam) );
%            mACF_Filter2(iLag,iParam) = corr( mParam_Accepted2(1:end-nMaxLag,iParam), mParam_Accepted2(iLag:end-nMaxLag+iLag-1,iParam) );
%        end
%     end
% end

mACF_Filter1 = zeros(nMaxLag+1,1);
mACF_Filter2 = zeros(nMaxLag+1,1);

for iParam = 1:nParamID 
    if vNotFixed(iParam,1) == 1
       mACF_Filter1 = [mACF_Filter1 zeros(nMaxLag+1,1)];
       mACF_Filter2 = [mACF_Filter2 zeros(nMaxLag+1,1)];
       for iLag = 1:nMaxLag+1    
           mACF_Filter1(iLag,end) = corr( mParam_Accepted1(1:end-nMaxLag,iParam), mParam_Accepted1(iLag:end-nMaxLag+iLag-1,iParam) );
           mACF_Filter2(iLag,end) = corr( mParam_Accepted2(1:end-nMaxLag,iParam), mParam_Accepted2(iLag:end-nMaxLag+iLag-1,iParam) );
       end
    end
end

mACF_Filter1 = mACF_Filter1(:,2:end);
mACF_Filter2 = mACF_Filter2(:,2:end);

%% Generate Scatter Plots - Part 2

% Panel 1: Lags 10 and 20

nIdLag1 = 11;
nIdLag2 = 21;

figure(1)
set(figure(1),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
hold on
scatter(mACF_Filter2(nIdLag1,:),mACF_Filter1(nIdLag1,:),100,'b','filled' )
scatter(mACF_Filter2(nIdLag2,:),mACF_Filter1(nIdLag2,:),150,'r','*' )
plot([0.5;1],[0.5;1],'Color','k','LineWidth',2)
%refline(1,0)
xlabel(strcat(nFilter2,''))
ylabel(strcat(nFilter1,''))
legend(['Lag ' num2str(nIdLag1-1)],['Lag ' num2str(nIdLag2-1)],'Location','northwest')
%legend('off')
%title(['Autocorrelations ' ])
box on
axis([0.5 1 0.5 1])
set(gca,'FontSize',30)
saveas(figure(1), [figPath '\ACF_Scatter_Lag_' num2str(nIdLag1-1) '_' num2str(nIdLag2-1) '.pdf' ] );
close all

% Panel 2: Lags 20 and 40

nIdLag1 = 31;
nIdLag2 = 41;

figure(1)
set(figure(1),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
hold on
scatter(mACF_Filter2(nIdLag1,:),mACF_Filter1(nIdLag1,:),100,'b','filled' )
scatter(mACF_Filter2(nIdLag2,:),mACF_Filter1(nIdLag2,:),150,'r','*' )
plot([0.5;1],[0.5;1],'Color','k','LineWidth',2)
%refline(1,0)
xlabel(strcat(nFilter2,''))
ylabel(strcat(nFilter1,''))
legend(['Lag ' num2str(nIdLag1-1)],['Lag ' num2str(nIdLag2-1)],'Location','northwest')
%legend('off')
%title(['Autocorrelations ' ])
box on
axis([0.5 1 0.5 1])
set(gca,'FontSize',30)
saveas(figure(1), [figPath '\ACF_Scatter_Lag_' num2str(nIdLag1-1) '_' num2str(nIdLag2-1) '.pdf' ] );
close all
