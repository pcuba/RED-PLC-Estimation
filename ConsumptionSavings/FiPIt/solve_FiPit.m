function Bdec_new = solve_FiPit(B_in, P_in, Z_in,Bdec,N_in)


[nz, nb] = size(Bdec);

STD_U = N_in(1);
RHO   = N_in(2);
GAMMAC= N_in(3);
R     = N_in(4);
M     = N_in(5);
BET  = N_in(6);

bmin = B_in(1);
bmax = B_in(nb);

[Zm, Bm ] = ndgrid(Z_in,B_in);

Cdec = Bdec - R*Bm + Zm;


for iz=1:nz
    for ib=1:nb
        
        b_use = B_in(ib);
        z_use = Z_in(iz);
        
        b_prime = Bdec(iz,ib);                     
        c_use = z_use - R*b_use + Bdec(iz,ib);
                
        if (b_prime < bmin)
            b_prime = bmin;
            c_use = b_prime - R*b_use + z_use;
        end
        


% allterp211 linear inter/extrapolation (2 states, 1 policies, 1 stoch comp)
% Inputs:
%   x*      :   Grid
%   x*i     :   Point to evaluate
%   pf*     :   Policy function
% Outputs:
%   o*      :   Interpolated/extrapolated values of dimension x*ipts
        
        % Bdec is in z x B space. 
        for iq=1:nz
            
            int_Ct_prime(iq) = allterp211(B_in, Z_in, Bdec(iq,ib),Z_in(iq), Cdec');
            
            %fastrap(nz,nb,Z_in,B_in,Z_in(iq),b_prime,Cd_in);
            
            int_Ct_prime(iq) = max(int_Ct_prime(iq),1E-20)^(-GAMMAC);
            
        end
                
        sol_Ct2 = P_in(iz,:)*(int_Ct_prime');
                
        Cdec_new(iz,ib) = (BET*R*sol_Ct2).^(-1/GAMMAC);
        
        Bdec_new(iz,ib) = Cdec_new(iz,ib) -z_use + R*b_use;
    
        if Bdec_new(iz,ib) > M*z_use
            Bdec_new(iz,ib) = M*z_use;
        end
        
    end
end
