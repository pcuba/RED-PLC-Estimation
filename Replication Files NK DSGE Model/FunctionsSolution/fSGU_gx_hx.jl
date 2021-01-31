
function fSGU_gx_hx(fy,fx,fyp,fxp);

#=
%[gx,hx,exitflag]=gx_hx(fy,fx,fyp,fxp,stake);
%computes the matrices gx and hx that define the first-order approximation to the solution
%of a dynamic stochastic general equilibrium model.
%Following the notation in Schmitt-Grohe and Uribe (JEDC, 2004), the model's equilibrium conditions
%take the form
%E_t[f(yp,y,xp,x)=0.
%The solution is of the form
%xp = h(x,sigma) + sigma * eta * ep
%y = g(x,sigma).
%The first-order approximations to the functions g and h around the point (x,sigma)=(xbar,0), where xbar=h(xbar,0), are:
%h(x,sigma) = xbar + hx (x-xbar)
%and
%g(x,sigma) = ybar + gx * (x-xbar),
%where ybar=g(xbar,0).
%The variable exitflag takes the values 0 (no solution), 1 (unique solution), 2 (indeterminacy), or 3 (z11 is not invertible).
%Inputs: fy fyp fx fxp stake
% The parameter stake ensures that all eigenvalues of hx are less than stake in modulus (the default is stake=1).
%Outputs: gx hx exitflag
%This program is a modified version of solab.m by Paul Klein  (JEDC, 2000).
%Modified on August 30 2009 by Ryan Chahrour (rc2374@columbia.edu) to replace qzdiv.m by ordqz.m  to increase speed.
%(c) Stephanie Schmitt-Grohe and Martin Uribe
%Date July 17, 2001, May 11 2006, November 19, 2009

=#

stake=1;
exitflag = 1;

NK = size(fx,2);

#Create system matrices A,B
A = [-fxp -fyp];
B = [fx fy];

#Complex Schur Decomposition
F = schur(A, B)
#s,t,q,z = schur(A, B)

#Pick non-explosive (stable) eigenvalues

slt = (abs.((F.beta)).<stake*abs.((F.alpha)));

nk=sum(slt);

#Reorder the system with stable eigs in upper-left
ordschur!(F,slt);

#Split up the results appropriately
z21 = F.Z[nk+1:end,1:nk];
z11 = F.Z[1:nk,1:nk];

s11 = F.S[1:nk,1:nk];
t11 = F.T[1:nk,1:nk];

#Identify cases with no/multiple solutions
if nk>NK
    println("The Equilibrium is Locally Indeterminate")
    exitflag=2;
elseif nk<NK
    println("No Local Equilibrium Exists")
    exitflag = 0;
end

if rank(z11)<nk;
    println("Invertibility condition violated")
    exitflag = 3;
end

if exitflag != 1
    
    return exitflag, nothing, nothing

else

    #Compute the Solution
    z11i = z11 \ Matrix{Float64}(I,nk,nk);
    gx = real(z21*z11i);
    hx = real(z11*(s11\t11)*z11i);

    return exitflag, gx,hx
end

end
