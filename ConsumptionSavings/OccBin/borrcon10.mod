
var b bnot c ec lb maxlev y ;
varexo eps_u ;

parameters RHO, BETA, M, R, STD_U, GAMMAC ;

model;
ec = c(1);
c = y + b - R*b(-1) ;
lb = 1/c^GAMMAC - BETA*R/c(+1)^GAMMAC ;
lb = 0;
bnot = M*y;
log(y) = RHO*log(y(-1)) + eps_u ;
maxlev = b-bnot;
end;



//shocks;
//var eps_u; stderr STD_U;
//end;












