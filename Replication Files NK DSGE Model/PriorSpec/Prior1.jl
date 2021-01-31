#--------------------------------------------------------------------------
# Script with prior information:
#
#   Each row has the following specification:
#
#      pdf,arg1,arg2,              mask,fix;
#      pdf shape of prior density
#                  1: BETA(mean,stdd)
#                  2: GAMMA(mean,stdd)
#                  3: NORMAL(mean,stdd)
#                  4: INVGAMMA(mode,df)
#                  5: UNIFORM(a,b)
#                  0: no prior
#      arg1    1st argument for pdf
#      arg2    2nd argument for pdf
#      Note:   for standard deviations, we start from
#              an Inverse Gamma distribution for variances,
#              where "mode" corresponds to the sqrt of the
#              mode of the distribution for the variance.
#              We then apply a change of variables to standard
#              deviations.
#      mask    1 if the parameter is fixed to constant
#      fix     fixed parameter value 
#
#--------------------------------------------------------------------------

mPriorSpec = [
2       2.00     0.20        0 2.0; # tau
1       0.10     0.05        0 0.2; # kappa
2       1.50     0.20        1 2.6; # psi1
2       0.75     0.20        1 0.98; # psi2
1       0.80     0.10        0 0.80; # rho_R
1       0.80     0.10        0 0.97; # rho_g
1       0.80     0.10        0 0.91; # rho_d
1       0.40     0.20        0 0.18; # rho_z
4       0.005    4.00        0 0.0016; # sigma_R
4       0.005    4.00        0 0.0025; # sigma_g
4       0.01     4.00        0 0.017; # sigma_d
4       0.01     4.00        0 0.0058; # sigma_z
2       1.00     0.50        1 0.72; # eta
2       0.20     0.10        1 0.10; # nu
2       1.00     0.10        1 1.00; # chi_h
2       1.20     0.20        0 1.2853; # gstar
2       1.00     0.40        0 0.0515; # rAnet
3       0.50     0.25        0 0.0190; # gamQnet
3       2.50     1.00        0 0.0363; # piAnet
3       2.50     1.00        1 0.5000; # pibarAnet #not used
3       0.00     0.002       0 0.0; # er0
3       0.00     0.012       0 0.0; # g0
3       0.00     0.008       0 0.0; # z0
3       0.00     0.17        0 0.0; # d0
3       0.00     0.02        1 0.0; # y0
3       0.00     0.03        0 0.0; # pi0
3       0.00     0.03        0 0.0; # c0
2       0.01     0.008       0 0.0; # R0star = R0 + R_SS
3       0.00     0.02        1 0.0] # ylag0
