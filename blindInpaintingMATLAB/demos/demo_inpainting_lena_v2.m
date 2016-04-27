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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (fractionmissing <= 0.70)&&(fractionmissing > 0.10)

    fprintf('Running Fast 2 phase TV...\n')
    t0 = cputime;
    [x_f2p, mask_f2p] = fastTV2Phase(y, 0.2, 1e-4);
    t_f2p = cputime-t0;

    PSNR_f2p = 10*log10(255^2*size(x_f2p,1)*size(x_f2p,2)/norm(x_f2p(:)-x0(:))^2);
    ISNR_f2p = 20*log10(norm(y(:)-x0(:))/norm(x_f2p(:)-x0(:)));
    SSIM_f2p = ssim(x0, x_f2p);
    numMaskErrors_f2p = sum(xor(mask(:),mask_f2p(:)));

    clc

    fprintf('\nFast 2 phase TV with %g per cent pixels missing and Gaussian noise sigma = %g .\n', fractionmissing*100, sigma)
    fprintf('%g seconds,\nISNR = %g dB, PSNR = %g dB, SSIM= %g,\nmask errors = %d, %% errors: %g %% \n',...
        t_f2p, ISNR_f2p, PSNR_f2p, SSIM_f2p, numMaskErrors_f2p, numMaskErrors_f2p*100/(M*N))

    figure, imagesc(x_f2p), colormap gray, axis image, axis off
    figure, imagesc(xor(mask, mask_f2p)), colormap gray, axis image, axis off

end      
            
close all
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (fractionmissing <= 0.70)

    fprintf('Running frame based inpainting...\n')

    t0 = cputime;
    [x_my, Dest, PSNR_my] = blindInpaintMingYan_v3(y/255, 640*5, 0, fractionmissing*100, 50, x0/255, 0);
    x_my = x_my*255;

    t_my = cputime - t0;

    ISNR_my = 20*log10(norm(y(:)-x0(:))/norm(x_my(:)-x0(:)));
    PSNR_my = 10*log10(255^2*size(x_my,1)*size(x_my,2)/norm(x_my(:)-x0(:))^2);
    SSIM_my = ssim(x0, x_my);
    numMaskErrors_my = sum(xor(mask(:),~Dest(:)));

    mask_my = ~Dest;

    fprintf('\nBlind inpainting (Ming Yan)\n')
    fprintf('%g seconds,\nISNR = %g dB, PSNR = %g dB, SSIM= %g,\nmask errors = %d, %% errors: %g %% \n',...
        t_my, ISNR_my, PSNR_my, SSIM_my, numMaskErrors_my, numMaskErrors_my*100/(M*N))

    figure, imagesc(x_my), colormap gray, axis image, axis off
    figure, imagesc(xor(mask, mask_my)), colormap gray, axis image, axis off
end

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (fractionmissing < 0.50)
    fprintf('Running kALS...\n')

    ALSopts.method = 2;
    ALSopts.iternum = 10;
    outPer = fractionmissing;

    t0 = cputime;
    [x_kals psnr_kals mask_kals D_kals] = KALS(x,y,10,1-outPer,0.01,[8 8],[1 1],ALSopts);
    x_kals(x_kals<0) = 0;
    x_kals(x_kals>255) = 255;
    t_kals = cputime - t0;

    PSNR_kals = 10*log10(255^2*size(x_kals,1)*size(x_kals,2)/norm(x_kals(:)-x0(:))^2);
    ISNR_kals = 20*log10(norm(y(:)-x0(:))/norm(x_kals(:)-x0(:)));
    SSIM_kals = ssim(x0, x_kals);
    numMaskErrors_kals = sum(xor(mask(:),mask_kals(:)));

    clc

    fprintf('\nkALS with %g per cent pixels missing and Gaussian noise sigma = %g .\n', fractionmissing*100, sigma)
    fprintf('%g seconds,\nISNR = %g dB, PSNR = %g dB, SSIM= %g,\nmask errors = %d, %% errors: %g %% \n',...
        t_kals, ISNR_kals, PSNR_kals, SSIM_kals, numMaskErrors_kals, numMaskErrors_kals*100/(M*N))

    figure, imagesc(x_kals), colormap gray, axis image, axis off
    figure, imagesc(xor(mask, mask_kals)), colormap gray, axis image, axis off

end

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bb=8; % block size
K=256; % number of atoms in the dictionary
N=256; 
iters = 1;
if fractionmissing >= 0.5
    iters = 1;
end

[x_ksvd, times_ksvd, mask_ksvd] = inpaintKSVD_v3(y_missing, sigma, N, K, bb, iters);

PSNR_ksvd = 10*log10(255^2*size(x_ksvd,1)*size(x_ksvd,2)/norm(x_ksvd(:)-x0(:))^2);
ISNR_ksvd = 20*log10(norm(y(:)-x0(:))/norm(x_ksvd(:)-x0(:)));
SSIM_ksvd = ssim(x0, x_ksvd);
numMaskErrors_ksvd = sum(xor(mask(:),mask_ksvd(:)));

fprintf('\KSVD with %g per cent pixels missing and Gaussian noise sigma = %g .\n', fractionmissing*100, sigma)
fprintf('%g seconds,\nISNR = %g dB, PSNR = %g dB, SSIM= %g,\nmask errors = %d, %% errors: %g %% \n',...
    times_ksvd(end), ISNR_ksvd, PSNR_ksvd, SSIM_ksvd, numMaskErrors_ksvd, numMaskErrors_ksvd*100/(M*N))

figure, imagesc(x_ksvd), colormap gray, axis image, axis off
figure, imagesc(xor(mask, mask_ksvd)), colormap gray, axis image, axis off

