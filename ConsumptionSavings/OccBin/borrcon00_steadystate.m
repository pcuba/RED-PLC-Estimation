
function [ys,check]=borrcon00_steadystate(junk,ys);

global M_

paramfile_borrcon00

nparams = size(M_.param_names,1);
for icount = 1:nparams
    eval(['M_.params(icount) = ',M_.param_names(icount,:),';'])
end


check=0;

b=M;
bnot=b;
maxlev=0;
c=1+M-R*M;
ec=c;
lb=(1-BETA*R)/c^GAMMAC;
y=1;


ys = [ b
bnot
c  
ec
lb
maxlev
y ] ;



