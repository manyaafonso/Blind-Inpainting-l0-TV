function objF = dataFidelity(x,y,model)

y = y(:);
x = x(:);

if model == 1
    objF = 0.5*norm(y-x)^2;
else
    if model == 2
        objF = sum( log(y+1e-4)-log(x+1e-4) - (y.^2)./(2*x) );
    else
        objF = sum( y.*log(x+1e-4)-x );
    end
end
    