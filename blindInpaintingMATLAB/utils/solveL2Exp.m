function x = solveL2Exp(r1, r2, beta, c, niters)

x = zeros(size(r1));

for iter = 1:niters
    x = x - ( x - r1 + beta*c*exp(c*x).*(exp(c*x)-r2) )./( 1 + beta*c*c*(2*exp(2*c*x) -c*r2.*exp(c*x) ) );
end
