function x = solveQuadratic1(a,b,c)

delta = sqrt(b.^2-4*a.*c);

x = min((-b+delta)./(2*a), (-b-delta)./(2*a) );

