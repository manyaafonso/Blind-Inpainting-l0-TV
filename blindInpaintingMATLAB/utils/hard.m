function y = hard(x,lambda)
y = double(abs(x)>lambda).*x;