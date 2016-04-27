function [x_est, mask_est, times] = blindInpaintAnscombe(y, tau1, tau2, maxiters, bias, tol )

[M,N] = size(y);

ya = Anscombe_forward(y);


ya = ya + bias;
g = log(ya);

t0 = cputime;
times(1) = 0;

u = zeros(M,N);
v = ones(M,N)*(-2);


% iters = 500;
% tau1 = 0.02; % 0.02 for 40% missing
% tau2 = 0.1;
% 
% u = zeros(size(g));
% %v = 01-mask)*(-10);
% v = ones(size(g))*(-10);

% figure;
x_hat = g;%exp(u);
x_hat_prev = x_hat;

mask_est = exp(v);
mask_est_prev = mask_est;

obj(1) = norm(ya(:)-x_hat(:).*mask_est(:))^2 + tau1*TVnorm(x_hat) + tau2*sum(v(:)~=0);


for iter = 1:maxiters
    
%     iter
    %u = soft(g-v,tau1);
    %u = chambolledenoise1d(g-v,tau1,100);
    %u = projk(g-v,tau1,5);
    u = tvdenoise(g-v,tau1,5);
    v = hard(g-u,tau2);   
    
    x_hat = exp(u);
    mask_est = exp(v);
    
    times(iter+1) = cputime-t0;
    
%     criterion1(iter) = norm(x_hat(:)-x_hat_prev(:))/norm(x_hat(:));
    %criterion2(iter) = norm(mask_est(:)-mask_est_prev(:))/norm(mask_est(:));
%     criterion2(iter) = sum(abs(mask(:)-mask_est(:))>=0.2);
    
    obj(iter+1) = norm(ya(:)-x_hat(:).*mask_est(:))^2 + tau1*TVnorm(x_hat) + tau2*sum(v(:)~=0);
    
    %if abs(obj(iter+1)-obj(iter))/obj(iter) <= 1e-10
    if ( norm(x_hat(:)-x_hat_prev(:))^2/norm(x_hat(:))^2 ) < tol
        break;
    end
    
    
    x_hat_prev = x_hat;
    mask_est_prev = mask_est;
    
    x_hat = exp(u);
%     imagesc(x_hat), colormap gray;
%     drawnow
    
%     tau1 = tau1/0.99;
%     tau2 = tau2/0.99;
end

% fprintf('%d iters, %g seconds.\n', iter, times(end))

x_hat = exp(u);
% figure, imagesc(exp(v)), colormap gray;
% % %a = double(g<=-6);
% figure, imagesc(x_hat), colormap gray;
% 
% figure, loglog(times, obj)
% figure, semilogy(times(2:end), criterion1)
% figure, semilogy(times(2:end), criterion2)

x_est = GenAnscombe_inverse_closed_form(x_hat,1e-4);
