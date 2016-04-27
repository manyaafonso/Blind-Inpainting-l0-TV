function [x_hat, mask_est, times] = denoiseRVImpulse_v2(y, tau1, tau2, maxiters, tol, bias, verbose)

y = y + bias;
[M,N] = size(y);

g = log10(y);


t0 = cputime;
times(1) = 0;

u = -4*ones(M,N);
%v = ones(size(mask))*(-4);
v = zeros(size(g));

figure;
x_hat = 10.^u;
x_hat_prev = x_hat;

mask_est = exp(v);
mask_est_prev = mask_est;

% obj(1) = norm(y(:)-x_hat(:).*mask_est(:))^2 + tau1*TVnorm(x_hat) + tau2*sum(v(:)~=0);

if verbose
    figure
end

for iter = 1:maxiters
    
    u = tvdenoisePG(g-v,tau1,5);
    v = soft(g-u,tau2);   
    
    x_hat = 10.^u;
    mask_est = double(v==0);
    
    times(iter+1) = cputime-t0;
    
    criterion = norm(x_hat(:)-x_hat_prev(:))/norm(x_hat(:));
    
    
    if ( criterion <= tol)
        break;
    end
%     
    
    x_hat_prev = x_hat;
    mask_est_prev = mask_est;
    
%     x_hat = exp(u);
    if verbose
            imagesc(x_hat), colormap gray;
        drawnow
    end
    
%     tau1 = tau1*1.02;
    %tau1 = min(tau1*1.1, 8);
%     tau1 = (M*N)/TVnorm(u);
%      tau2 = min(tau2*1.02, 0.1);

    tau1 = min(tau1/0.99, 1);
    
    tau2 = max(tau2*0.999, 0.01);
    
    
end


 