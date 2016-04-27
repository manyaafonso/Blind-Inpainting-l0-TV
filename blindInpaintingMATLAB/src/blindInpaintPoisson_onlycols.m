function [x_hat, mask_est, times, err_x, err_mask] = blindInpaintPoisson_onlycols(y, lambda1, lambda2, mu1, mu2,...
    maxiters, tol, chambolleit, verbose, x_true, mask_true)

if nargin < 11
    compute_err_mask = 0;
else
    compute_err_mask = 1;
end

if nargin < 10
   compute_err_x = 0;
else
   compute_err_x = 1;
end 

if nargin < 9
    verbose = 0;
end

[M,N] = size(y);

t0 = cputime;
times(1) = 0;

u = zeros(M,N);
v = ones(M,N)*(-4);

x_hat = y;
mask_est = double(exp(v)>0.2);

z = 0*y;
w = 0*y;
a = 0*ones(size(y));

dz = 0*z;
dw = 0*w;
% df = 0*f;
if compute_err_x
	err_x(1) = norm(x_hat(:)-x_true(:))^2/(M*N);
end

if compute_err_mask
	err_mask(1) = sum(xor(mask_est(:),mask_true(:)));
end

if verbose
figure
end
for iter = 1:maxiters
    if verbose
        iter
    end
    
    %z = tvdenoisePG(u-dz, lambda1/mu1, chambolleit);
    z = zeros(M,N);
    for i = 1:N
        z(:,i) = projk1d(u(:,i)-dz(:,i), lambda1/mu1, chambolleit);
    end
        
    
    w = hard(v-dw, lambda2/mu2);
    
    u = solveLinPlusExp(mu1*ones(M,N), exp(v), ones(M,N), -mu1*(z+dz)-y );
    
    v = solveLinPlusExp(mu2*ones(M,N), exp(u), ones(M,N), -mu2*(w+dw)-y );
    

    x_hat_prev = x_hat;

    x_hat = exp(u);
    mask_est = double(exp(v)>0.5);
    
    if ( norm(x_hat(:)-x_hat_prev(:))^2/norm(x_hat(:))^2 ) < tol
        break
    end
    
    dz = dz - u + z;
    dw = dw - v + w;
    
    times(iter+1) = cputime - t0;
    
    if compute_err_x
        err_x(iter+1) = norm(x_hat(:)-x_true(:))^2/(M*N);
    end
    
    if compute_err_mask
    	err_mask(iter+1) = sum(xor(mask_est(:),mask_true(:)));
    end

    if verbose
    imagesc(x_hat), colormap gray;
    drawnow
    end
end

