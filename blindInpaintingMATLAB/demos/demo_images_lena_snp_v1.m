close all;
clear;

%addpath(genpath('/mnt/home/biomedic/work/inpaintingAndInpulseNoise/blindInpaintingMATLAB/'))
addpath(genpath('..'))

example = 'Lena';

x0 = imresize(imread('lena.png'), 0.5);
[M,N] = size(x0);
    
x0 = 255*im2double(x0); 
x = x0;

Psig  = norm(x0(:))^2/(M*N);

gSNR = 10;
sigma = norm(x0(:)-mean(x0(:)))/sqrt(N*M*10^(gSNR/10));
   
        
fractionmissing = 0.5;

yg = x0 + sigma*randn(size(x0));
yg = max(yg,0);

mask = double(rand(M,N)>fractionmissing);
snpnoise0 = 255*double(rand(M,N)>0.5);

y = yg.*mask + (1-mask).*snpnoise0;

figure, imagesc(y), colormap gray, axis image, axis off

%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('S&P denoising with %g per cent pixels missing and sigma = %g.\n',fractionmissing*100, sigma)
%fprintf('%d iterations out of %d.\n', mc, MCiters)
fprintf('Running Proposed method...\n')

bias = 1e-3;
maxiters = 1000;
tau1 = 1; % 0.005
tau2 = 0.1;%0.05

[x_hat, mask_est, times] = denoiseSnP_v1(y, tau1, tau2, maxiters, 1e-4, bias, 0);

ISNR = 20*log10(norm(y(:)-x0(:))/norm(x_hat(:)-x0(:)));
PSNR = 10*log10(255^2/(norm(x_hat(:)-x0(:))^2/numel(x)));
SSIM = ssim(x0, x_hat);

figure, imagesc(x_hat), colormap gray, axis image, axis off

fprintf('%d iters, %g seconds, ISNR = %g dB, PSNR = %g dB, SSIM= %g \n', length(times)-1, times(end), ISNR, PSNR, SSIM)

            