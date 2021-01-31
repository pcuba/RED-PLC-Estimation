function [depvar, ylab] = fTransform(DR,vname,par,yunits);
    
    eval(['depvar = DR.'  vname ';']);
    eval(['varss  = par.' vname 'ss;']);
    
    if yunits == 0;
        % No transformation
        depvar = depvar;
        ylab = 'Level';
    elseif yunits == 1;
        % Deviation from SS
        depvar = log(depvar/varss)*100;
        ylab = '$\% \Delta $ from S.S.';
    elseif yunits ==2;
        % Annualized rate of gross variable
        depvar = 400*log(depvar);
        ylab = '$\% $ annualized';
    elseif yunits ==3;
        % Annualized rate of rate variable
        depvar = 400*(depvar);
        ylab = '$\% $ annualized'; 
    end
    
end

