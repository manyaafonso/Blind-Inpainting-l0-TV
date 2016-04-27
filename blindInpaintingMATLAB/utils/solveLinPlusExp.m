function x = solveLinPlusExp(a,b,c,d)

numIters = 5;

x = zeros( size(a) );

for iter = 1:numIters
    
    x = x - ( a.*x + b.*exp(c.*x) + d )./( a + b.*c.*exp( c.*x ) );
end