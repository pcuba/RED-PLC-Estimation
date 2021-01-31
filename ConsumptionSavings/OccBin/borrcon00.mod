
var b bnot c ec lb maxlev y ;
varexo eps_u ;

parameters RHO, BETA, M, R, STD_U, GAMMAC  ;

R    = 1.05;
BETA = 0.945;
RHO  = 0.9;
STD_U= 0.010;
M = 1;
GAMMAC = 1; 

model;
ec = c(1);
c = y + b - R*b(-1) ;
b =  M*y ;
bnot =  M*y ;
lb  = 1/c^GAMMAC - BETA*R/c(+1)^GAMMAC ;
log(y) = RHO*log(y(-1)) + eps_u ;
maxlev = b-bnot;
end;





//shocks;
//var eps_u; stderr STD_U;
//end;

//steady;




stoch_simul(order=1,noprint,nomoments,irf=0);









