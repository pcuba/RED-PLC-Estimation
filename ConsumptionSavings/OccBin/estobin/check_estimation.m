function mmm = check_estimation(foldername,modfilename,cutoff1)

% Example: 
% check_estimation({'C:\Dropbox\E\occbin_Estimation\borrcon','C:\Dropbox\E\occbin_Estimation\borrcon_alt1'},'borrcon00',15)

do_two=0;
set(0,'DefaultLineLineWidth',2)
warning off
format short

nx=size(foldername,2); 

for ix=1:nx
    folder=foldername{ix}; 
    
    disp(' ')
    disp(' ')
    disp(folder)
    
    eval([ 'load ' folder  '\datavec '])
    eval([ 'load ' folder  '\' modfilename '_results M_ '])
    eval([ 'load ' folder  '\params_matrix params_matrix '])
    
    poststory=-fstory(:,1);  
    [ maxl imaxl ] = max(poststory) ;


    xstory_copy=xstory;
    poststory_copy=-fstory(:,1);
    xstory(poststory<max(poststory)-cutoff1,:) = NaN ;
    poststory(poststory<max(poststory)-cutoff1)=NaN;
    
    disp('PARAMETER   MAXIMUM    INITIAL')
    npar=numel(params_labels);
    for ii=1:npar
        trspaces=blanks(10-size(char(params_labels(ii,:)),2));
        trspaces2=blanks(11-size(char(params_matrix(ii,6)),2));
        disp([ char(params_labels(ii,:))  trspaces ' = ' num2str(xstory(imaxl,ii),'%0.4f') ';    '  num2str(xstory_copy(1,ii),'%0.4f') ])
    end
    disp([ 'Posterior = ' num2str(max(poststory),8)  ';    '  num2str(poststory_copy(1),8)])
    
    
    
    figure
    if npar<10; nr=3; nc=3; end
    if npar>10; nr=3; nc=7; end
    if npar>21; nr=4; nc=6; end
    for ii=1:npar
        subplot(nr,nc,ii)
        plot(xstory(:,ii),poststory,'.'); hold on
        plot(xstory(1,ii),poststory(1),'or'); hold on
        plot(xstory(imaxl,ii),poststory(imaxl),'og'); hold on
        %         ylim([poststory(1) max(poststory) ])
        axis tight
        title(params_labels(ii,:))
    end
    
    subplot(nr,nc,ii+1)
    plot(1:numel(poststory),poststory,'.'); hold on
    plot(1,poststory(1),'or'); hold on
    plot(imaxl,poststory(imaxl),'og'); hold on
    %     ylim([poststory(1) max(poststory) ])
    if numel(poststory)>1
    xlim([1 numel(poststory)])
    end
    title('Posterior (or likelihood)')
    
    disp('Mode stored in aaa')
    aaa(1:npar,ix)=xstory(imaxl,:)    ;
    bbb(ix)=poststory(imaxl);
    ccc(ix)=length(datavec);
    ddd(ix)=imaxl;
    
    
    if size(xstory,2)==2 || do_two==1
        
        if do_two==1
            disp('find two variables that move')
            ivars=find(nanmean(abs(diff(xstory)))>0);
            ivar1=ivars(1);
            ivar2=ivars(2);
        else
            ivar1=1;
            ivar2=2;
        end
            
        
        [z]=poststory(poststory>max(poststory)-cutoff1);
        [xcopy]=xstory(find(poststory>max(poststory)-cutoff1),:);
        
        % Data
        xx=xcopy(:,ivar1);
        y=xcopy(:,ivar2);
        xx(isnan(xx))=[];
        y(isnan(y))=[];
        z(isnan(z))=[];
        
        % Interpolated points
        nn = 31;
        xi = linspace(min(xx), max(xx), nn);
        yi = linspace(min(y), max(y), nn);
        [xi, yi] = meshgrid(xi, yi);
        
        % Plots
        figure
        
        subplot(3,1,1)
        frequency=1+round(500*(max(poststory)-poststory)/max(poststory));
        scatter(xstory(:,ivar1),xstory(:,ivar2),frequency,'filled'); hold on
        plot(xstory(imaxl,ivar1),xstory(imaxl,ivar2),'*','color','r');        
        for i=1:numel(xstory(:,ivar1))
            text(xstory(i,ivar1),xstory(i,ivar2),num2str(max(poststory)-poststory(i),'%0.0f')); hold on
        end


        subplot(3,1,2)
        interpolator=TriScatteredInterp(xx,y,z,'linear');
        zi=interpolator(xi,yi);
        mesh(xi, yi, zi)
        hold on
        plot3(xx,y,z,'.','MarkerSize',15) %nonuniform
        xlim([min(xx) max(xx)])
        ylim([min(xx) max(xx)])
        xlabel(params_labels(ivar1))
        ylabel(params_labels(ivar2))
        axis tight

        

        subplot(3,1,3)
        [c,h] = contour(xi,yi,(zi));
        set(h,'ShowText','on','TextStep',get(h,'LevelStep')*.25)
        colormap cool
        xlabel(params_labels(ivar1))
        ylabel(params_labels(ivar2))

        
        
    end
    
end






if ix==2
    disp(' ')
    disp('Comparison')
    for ii=1:npar
        trspaces=blanks(10-size(char(params_labels(ii,:)),2));
        trspaces1=blanks(3);
        disp([ char(params_labels(ii,:)) ...
            trspaces ' = ' num2str(aaa(ii,1),'%0.4f')  ...
            trspaces1   num2str(aaa(ii,2),'%0.4f')   ])
    end
    disp([ 'fval         ' num2str(bbb(1),'%0.8f')  ...
        trspaces1   num2str(bbb(2),'%0.8f') ])
    disp([ 'iters         ' num2str(ccc(1),'%0.0f')  ...
        trspaces1 trspaces1    num2str(ccc(2),'%0.0f') ])
    disp([ 'max @         ' num2str(ddd(1),'%0.0f')  ...
        trspaces1 trspaces1    num2str(ddd(2),'%0.0f') ])
end
if ix==3
    disp(' ')
    disp('Comparison')
    for ii=1:npar
        trspaces=blanks(10-size(char(params_labels(ii,:)),2));
        trspaces1=blanks(3);
        disp([ char(params_labels(ii,:)) ...
            trspaces ' = ' num2str(aaa(ii,1),'%0.4f')  ...
            trspaces1   num2str(aaa(ii,2),'%0.4f')  ...
            trspaces1   num2str(aaa(ii,3),'%0.4f') ])
    end
    disp(' ')
    disp([ 'fval         ' num2str(bbb(1),'%0.1f')  ...
        trspaces1   num2str(bbb(2),'%0.1f')  ...
        trspaces1   num2str(bbb(3),'%0.1f') ])
    disp([ 'iters         ' num2str(ccc(1),'%0.0f')  ...
        trspaces1 trspaces1    num2str(ccc(2),'%0.0f') ...
        trspaces1 trspaces1    num2str(ccc(3),'%0.0f') ]);
    disp([ 'max @         ' num2str(ddd(1),'%0.0f')  ...
        trspaces1 trspaces1    num2str(ddd(2),'%0.0f') ...
        trspaces1 trspaces1    num2str(ddd(3),'%0.0f') ]);
end
if ix==4
    disp(' ')
    disp('Comparison')
    for ii=1:npar
        trspaces=blanks(10-size(char(params_labels(ii,:)),2));
        trspaces1=blanks(3);
        disp([ char(params_labels(ii,:)) ...
            trspaces ' = ' num2str(aaa(ii,1),'%0.4f')  ...
            trspaces1   num2str(aaa(ii,2),'%0.4f')  ...
            trspaces1   num2str(aaa(ii,3),'%0.4f')  ...
            trspaces1   num2str(aaa(ii,4),'%0.4f') ])
    end
    disp(' ')
    disp([ 'fval       ' num2str(bbb(1),'%0.2f')  ...
        trspaces1   num2str(bbb(2),'%0.2f')  ...
        trspaces1   num2str(bbb(3),'%0.2f')  ...
        trspaces1   num2str(bbb(4),'%0.2f') ])
    disp([ 'iters         ' num2str(ccc(1),'%0.0f')  ...
        trspaces1 trspaces1    num2str(ccc(2),'%0.0f') ...
        trspaces1 trspaces1    num2str(ccc(3),'%0.0f') ...
        trspaces1 trspaces1    num2str(ccc(4),'%0.0f') ]);
    disp([ 'max @         ' num2str(ddd(1),'%0.0f')  ...
        trspaces1 trspaces1    num2str(ddd(2),'%0.0f') ...
        trspaces1 trspaces1    num2str(ddd(3),'%0.0f') ...
        trspaces1 trspaces1    num2str(ddd(4),'%0.0f') ]);
end
if ix==5
    disp(' ')
    disp('Comparison')
    for ii=1:npar
        trspaces=blanks(10-size(char(params_labels(ii,:)),2));
        trspaces1=blanks(3);
        disp([ char(params_labels(ii,:)) ...
            trspaces ' = ' num2str(aaa(ii,1),'%0.4f')  ...
            trspaces1   num2str(aaa(ii,2),'%0.4f')  ...
            trspaces1   num2str(aaa(ii,3),'%0.4f')  ...
            trspaces1   num2str(aaa(ii,4),'%0.4f')  ...
            trspaces1   num2str(aaa(ii,5),'%0.4f') ])
    end
    disp(' ')
    disp([ 'fval       ' num2str(bbb(1),'%0.2f')  ...
        trspaces1   num2str(bbb(2),'%0.2f')  ...
        trspaces1   num2str(bbb(3),'%0.2f')  ...
        trspaces1   num2str(bbb(4),'%0.2f')  ...
        trspaces1   num2str(bbb(5),'%0.2f') ])
    disp([ 'iters         ' num2str(ccc(1),'%0.0f')  ...
        trspaces1 trspaces1    num2str(ccc(2),'%0.0f') ...
        trspaces1 trspaces1    num2str(ccc(3),'%0.0f') ...
        trspaces1 trspaces1    num2str(ccc(4),'%0.0f') ...
        trspaces1 trspaces1    num2str(ccc(5),'%0.0f') ]);
    disp([ 'max @         ' num2str(ddd(1),'%0.0f')  ...
        trspaces1 trspaces1    num2str(ddd(2),'%0.0f') ...
        trspaces1 trspaces1    num2str(ddd(3),'%0.0f') ...
        trspaces1 trspaces1    num2str(ddd(4),'%0.0f') ...
        trspaces1 trspaces1    num2str(ddd(5),'%0.0f') ]);
end
if ix==6
    disp(' ')
    disp('Comparison')
    for ii=1:npar
        trspaces=blanks(10-size(char(params_labels(ii,:)),2));
        trspaces1=blanks(3);
        disp([ char(params_labels(ii,:)) ...
            trspaces ' = ' num2str(aaa(ii,1),'%0.4f')  ...
            trspaces1   num2str(aaa(ii,2),'%0.4f')  ...
            trspaces1   num2str(aaa(ii,3),'%0.4f')  ...
            trspaces1   num2str(aaa(ii,4),'%0.4f')  ...
            trspaces1   num2str(aaa(ii,5),'%0.4f')  ...
            trspaces1   num2str(aaa(ii,6),'%0.4f') ])
    end
end



% CRAZY ATTEMPT TO GET HESSIAN
% my y is likstory
% my x is xstory
% Regress y on constant, x and x-squared


% [ hessian_reg stdh_reg hessian_fmin stdh_fmin ] = compute_hessian(xstory_copy,fstory,10);
%  [ stdh_reg./ stdh_fmin ]

% if ix==1
%     disp('PARAMETER   MODE     STD_REG    STD_FMIN   RATIO')
%     for ii=1:npar
%         trspaces=blanks(10-size(char(params_labels(ii,:)),2));
%         trspaces1=blanks(3);
%         disp([ char(params_labels(ii,:)) ...
%             trspaces ' = ' num2str(aaa(ii,1),'%0.4f')  ...
%             trspaces1   num2str(stdh_reg(ii,1),'%0.4f') ...
%             '; ' trspaces1   num2str(stdh_fmin(ii,1),'%0.4f') ...
%             '; ' trspaces1   num2str(stdh_fmin(ii,1)/stdh_reg(ii,1),'%0.4f')  ])
%     end
% end










% for i=1:6; subplot(3,2,i); plot(sequence(:,i)); hold on; plot(filtered_errs(:,i),'r'); end