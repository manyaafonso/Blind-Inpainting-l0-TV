function [x_hat, mask_est, obj, times, err_x, err_mask] = blindInpaint_impulsenoise_v1(varargin)

if nargin < 3
    error('Must specify observed image y, and regularization parameters tau1, tau2!')
else
    y = varargin{1};
    tau1 = varargin{2};
    tau2 = varargin{3};
        
    % default parameters
    bias = 1e-3;
    compute_errx = 0;
    compute_errm = 0;
    maxiters = 500;
    chambolleit = 5;
    verbose = 0;
    stopcriterion = 0;
    cont_tau1 = 0;
    cont_tau2 = 0;

    if nargin == 4
        %%% if options are provided, use them, else use default parameters.
        opts = varargin{4};
        if isfield(opts,'bias')
            bias = opts.bias;
        end
        if isfield(opts,'maxiters')
            maxiters = opts.maxiters;
        end
        if isfield(opts,'chambolleit')
            chambolleit = opts.chambolleit;
        end
        if isfield(opts,'x_true')
            compute_errx = 1;
            x_true = opts.x_true;
        end
        if isfield(opts,'mask_true')
            compute_errm = 1;
            mask_true = opts.mask_true;
        end
        if isfield(opts,'verbose')
            verbose = opts.verbose;
        end
        if isfield(opts,'stopcriterion')
            stopcriterion = opts.stopcriterion;
            if isfield(opts,'tol')
                tol = opts.tol;
            else
                if stopcriterion
                    error('Must specify tolerance for stopping criterion.');
                end
            end
        end
        if isfield(opts,'continue_tau1')
            cont_tau1 = 1;
            updateRule_tau1 = opts.continue_tau1;
        end
        if isfield(opts,'continue_tau2')
            cont_tau2 = 1;
            updateRule_tau2 = opts.continue_tau2;
        end

    end
end

[M,N] = size(y);

%%% add positive bias term to avoid zeros
y = y + bias;

g = log(y);

t0 = cputime;
times(1) = 0;

u = zeros(size(g));
v = ones(M,N)*(log10(bias)-1);

x_hat = exp(u);
x_hat_prev = x_hat;

mask_est = exp(v);


if compute_errx
    err_x(1) = norm(x_true(:)-x_hat(:))^2/(M*N);
end

if compute_errm
    err_mask(1) = norm(mask_true(:)-mask_est(:))^2/(M*N);
end

%obj(1) = norm(y(:)-x_hat(:).*mask_est(:))^2 + tau1*TVnorm(x_hat) + tau2*sum(v(:)~=0);
%obj(1) = norm(y(:)-x_hat(:).*mask_est(:))^2 + tau1*TVnorm(x_hat) + tau2*sum(mask_est(:)~=0);
obj(1) = norm(g(:)-u(:)-v(:))^2 + tau1*TVnorm(u) + tau2*sum(abs(v(:)));

if verbose
    figure
end
for iter = 1:maxiters
    
    u = projk(g-v,tau1,chambolleit);
    v = soft(g-u,tau2);   
    
    x_hat = exp(u);
    mask_est = exp(v);
    
    if verbose
        imagesc(x_hat), colormap gray;
        drawnow
    end
    
    times(iter+1) = cputime-t0;
    
    %obj(iter+1) = norm(y(:)-x_hat(:).*mask_est(:))^2 + tau1*TVnorm(x_hat) + tau2*sum(v(:)~=0);
    %obj(iter+1) = norm(y(:)-x_hat(:).*mask_est(:))^2 + tau1*TVnorm(x_hat) + tau2*sum(mask_est(:)~=0);
    obj(iter+1) = norm(g(:)-u(:)-v(:))^2 + tau1*TVnorm(u) + tau2*sum(abs(v(:)));
    
    
    if compute_errx
        err_x(iter+1) = norm(x_true(:)-x_hat(:))^2/(M*N);
    end

    if compute_errm
        err_mask(iter+1) = norm(mask_true(:)-mask_est(:),1);
    end

    
    if stopcriterion
        
        switch stopcriterion
            case 1
                criterion = abs( (obj(iter+1)-obj(iter))/obj(iter) );
            case 2
                criterion = norm(x_hat(:)-x_hat_prev(:))/norm(x_hat(:));
            case 3
                criterion = obj(t+1);
            otherwise
                error('Invalid stopping criterion!')
        end
        
        if criterion < tol
            if verbose
                fprintf('Convergence reached.\n')
            end
            break;
        end
        
    end
    
    x_hat_prev = x_hat;
    
    if cont_tau1
        tau1 = updateRule_tau1(tau1);
    end
    
    if cont_tau2
        tau2 = updateRule_tau2(tau2);
    end
    
    
end




