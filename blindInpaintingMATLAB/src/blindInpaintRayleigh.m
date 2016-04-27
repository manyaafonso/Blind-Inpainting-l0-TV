function [x_hat, mask_est, times, err_x, err_mask] = blindInpaintRayleigh(y, lambda1, lambda2, mu1, mu2,...
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

g = log(y);

u = zeros(M,N);
v = 0*ones(M,N)*(-4);

x_hat = exp(u);
x_hat_prev = x_hat;


mask_est = exp(v);
mask_est_prev = mask_est;

obj(1) = norm(y(:)-x_hat(:).*mask_est(:))^2 + lambda1*TVnorm(x_hat) + lambda2*sum(v(:)~=0);

ysq = y.^2;

dz = zeros(M,N);
dw = zeros(M,N);

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
    %z = projk(u-dz,tau1/mu_u,20);
    z = tvdenoisePG(u-dz,lambda1/mu1,chambolleit);
    u = solveLinPlusExp( mu1*ones(M,N), -0.25*ysq.*exp(-v), -ones(M,N), ones(M,N)-mu1*(z+dz) );
    
	w = hard(v-dw,lambda2/mu2);
    
    v = solveLinPlusExp( mu2*ones(M,N), -0.25*ysq.*exp(-u), -ones(M,N), ones(M,N)-mu2*(w+dw) );
    
    
    
    dz = dz - (u-z);
    dw = dw - (v-w);
    
    
    x_hat = exp(u);
    mask_est = exp(v);
    
    times(iter+1) = cputime-t0;
    
    if ( norm(x_hat(:)-x_hat_prev(:))^2/norm(x_hat(:))^2 ) < tol
        break;
    end
    
    
    x_hat_prev = x_hat;
    mask_est_prev = mask_est;
    
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


x_hat = exp(u);
