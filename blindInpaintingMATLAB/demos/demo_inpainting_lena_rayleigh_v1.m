close all
clear all

%addpath(genpath('/mnt/home/biomedic/work/inpaintingAndInpulseNoise/blindInpaintingMATLAB/'))
addpath(genpath('..'))

example = 'Lena';

x0 = double(rgb2gray(imread('lena512color.tiff')));
[M,N] = size(x0);

x = abs(x0)/max(x0(:));

[M,N] = size(x0);
    


fractionmissing = 0.5;


noise0 = sqrt(raylrnd(ones(M,N)));
y0 = noise0.*x;

mask = double(rand(M,N)>fractionmissing);

y = y0.*mask;

figure, imagesc(y), colormap gray, axis image, axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Blind inpainting with Rayleigh noise and %g per cent pixels missing.\n', fractionmissing*100)
fprintf('Running Proposed method...\n')

bias = 1e-3;
iters = 2500;
tol = 1e-6;
chambolleit = 10;

tau1 = 5; % 0.005
tau2 = 1;
mu1 = 1;
mu2 = 1;

[x_hat, mask_est, times, err_x, err_mask] = blindInpaintRayleigh(y+bias, tau1, tau2, mu1, mu2,...
    iters, tol, chambolleit, 0, x, mask);
fprintf('%d iters, %g seconds.\n', length(times)-1, times(end))

MAE = sum(abs(x(:)-x_hat(:)))/numel(x);
MAEn = MAE/sum(abs(x(:)));
SSIM = ssim(x, x_hat);
numMaskErrors = sum(xor(mask(:),double(mask_est(:)>0.1)));

fprintf('Blind inpainting with Rayleigh noise and %g per cent pixels missing.\n',fractionmissing*100)
fprintf('%d iters, %g seconds, MAE = %g , MAEn = %g , SSIM= %g, mask errors = %g pc \n',...
    length(times)-1, times(end), MAE, MAEn, SSIM, 100*numMaskErrors/(M*N))

figure, imagesc(x_hat), colormap gray, axis image, axis off
figure, imagesc(xor(mask,double(mask_est>0.5))), colormap gray, axis image, axis off

            