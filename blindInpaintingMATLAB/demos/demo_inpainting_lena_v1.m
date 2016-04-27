close all
clear all

%addpath(genpath('/mnt/home/biomedic/work/inpaintingAndInpulseNoise/blindInpaintingMATLAB/'))
addpath(genpath('..'))

x0 = imresize(imread('lena.png'), 0.5);
[M,N] = size(x0);

x0 = 255*im2double(x0); 
x = x0;

gSNR = 10;
fractionmissing = 0.75;


Psig  = norm(x0(:))^2/(M*N);
sigma = norm(x0(:)-mean(x0(:)))/sqrt(N*M*10^(gSNR/10));
yg = x0 + sigma*randn(size(x0));
yg = max(yg,0);


mask = double(rand(M,N)>fractionmissing);

y = yg.*mask;

figure, imagesc(y), colormap gray, axis image, axis off

fprintf('\nLena image, blind inpainting with Gaussian noise SNR = %g dB and %g per cent pixels missing.\n',...
    gSNR, fractionmissing*100)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
bias = 1e-3;
updateTau1 = @(x) min(x*1.02, 2);
updateTau2 = @(x) min(x*1.02, 4);
maxiters = 500;
tau1 = 0.008; % 0.005
tau2 = 0.1;
y_missing = max(y,0);


opts = struct('maxiters',500,'chambolleit',5,'stopcriterion',2,...
    'tol', 1e-3,'verbose',0,'x_true',x,'mask_true',mask,...
    'continue_tau1',updateTau1, 'continue_tau2',updateTau2);
[x_hat, mask_est, obj, times, err_x, err_mask] = blindInpaint_v2(y_missing, tau1, tau2, opts);
x_hat = real(x_hat);
fprintf('%d iters, %g seconds.\n', length(times)-1, times(end))

ISNR = 20*log10(norm(y(:)-x0(:))/norm(x_hat(:)-x0(:)));
PSNR = 10*log10(255^2/err_x(end));
SSIM = ssim(x0, x_hat);
numMaskErrors = sum(xor(mask(:),double(mask_est(:)>0.5)));

fprintf('\nBlind inpainting with Gaussian noise SNR = %g dB and %g per cent pixels missing.\n',gSNR, fractionmissing*100)
fprintf('%d iters, %g seconds,\nISNR = %g dB, PSNR = %g dB, SSIM= %g,\nmask errors = %d, %% errors: %g %% \n',...
    length(times)-1, times(end), ISNR, PSNR, SSIM, numMaskErrors, numMaskErrors*100/(M*N))

figure, imagesc(x_hat), colormap gray, axis image, axis off

figure, imagesc(xor(mask,double(mask_est>0.5))), colormap gray, axis image, axis off

