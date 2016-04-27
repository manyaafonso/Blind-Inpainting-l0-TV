close all
clear all

%addpath(genpath('/mnt/home/biomedic/work/inpaintingAndInpulseNoise/blindInpaintingMATLAB/'))
addpath(genpath('..'))

example = 'Lena';

x0 = 0.5*double(rgb2gray(imread('lena512color.tiff')));
[M,N] = size(x0);
    
%y0 = poissrnd(x0);

fractionmissing = 0.5;

y0 = poissrnd(x0);

mask = double(rand(M,N)>fractionmissing);

y = y0.*mask;

figure, imagesc(y), colormap gray, axis image, axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Blind inpainting with Poisson noise and %g per cent pixels missing.\n', fractionmissing*100)
fprintf('Running Proposed method...\n')

lambda1 = 500; % 0.005 200 for 0.5, 500 for 0.75
lambda2 = 200;
mu1 = 50;
mu2 = 100;

[x_hat, mask_est, times, err_x, err_mask] = blindInpaintPoisson(y, lambda1, lambda2, mu1, mu2,...
        500, 1e-6, 10, 0, x0, mask);
fprintf('%d iters, %g seconds.\n', length(times)-1, times(end))

MAE = sum(abs(x0(:)-x_hat(:)))/numel(x0);
MAEn = MAE/sum(abs(x0(:)));
SSIM = ssim(x0, x_hat);
numMaskErrors = sum(xor(mask(:),double(mask_est(:)>0.1)));


fprintf('Blind inpainting with Poisson noise and %g per cent pixels missing.\n',fractionmissing*100)
fprintf('%d iters, %g seconds, MAE = %g , MAEn = %g , SSIM= %g, mask errors = %g pc \n',...
    length(times)-1, times(end), MAE, MAEn, SSIM, 100*numMaskErrors/(M*N))

figure, imagesc(x_hat), colormap gray, axis image, axis off
figure, imagesc(xor(mask,double(mask_est>0.5))), colormap gray, axis image, axis off

