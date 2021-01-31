function o1 = allterp211(  x1, x2, ...
                           x1i, x2i, ...
                           pf1)

% allterp211 linear inter/extrapolation (2 states, 1 policies, 1 stoch comp)
% Inputs:
%   x*      :   Grid
%   x*i     :   Point to evaluate
%   pf*     :   Policy function
% Outputs:
%   o*      :   Interpolated/extrapolated values of dimension x*ipts

% Grid lengths
nx1 = length(x1);
nx2 = length(x2);

% Number of stochastic realizations
x2ipts = length(x2i);

% Preallocate output
o1 = zeros(x2ipts,1);

s1 = x1(2) - x1(1);
x1i_min = x1i - x1(1);
loc1 = min(nx1-1,max(1,floor(x1i_min/s1) + 1));

for i2 = 1:x2ipts
    s2 = x2(2) - x2(1);
    x2i_min = x2i(i2) - x2(1);
    loc2 = min(nx2-1,max(1,floor(x2i_min/s2) + 1));

    xi = [x1i x2i(i2)];
    xi_left = [x1(loc1) x2(loc2)];
    xi_right = [x1(loc1+1) x2(loc2+1)];

    w_2 = (xi - xi_left)./ (xi_right - xi_left);
    w_1 = 1 - w_2;
    w1 = [w_1(1) w_2(1)];
    w2 = [w_1(2) w_2(2)];
    
    for m2 = 0:1
        for m1 = 0:1
            o1(i2) = o1(i2) + w1(m1+1)*w2(m2+1)*pf1(loc1+m1,loc2+m2);
        end
    end
end