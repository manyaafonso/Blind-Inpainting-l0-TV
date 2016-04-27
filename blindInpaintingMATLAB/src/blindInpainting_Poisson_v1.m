close all;
clear;

% addpath /home/manya/work/SALSA_v2.0/utils
% addpath(genpath('../../'))
%addpath(genpath('../../../blindInpaintingMATLAB'))
addpath(genpath('../../src'))
addpath(genpath('../../imgs'))
addpath(genpath('../../utils'))


%x0 = double(imread('cameraman.tif'));
 x0 = double(rgb2gray(imread('lena512color.tiff')));
x0 = 0.5*imresize(x0,0.5);

[M,N] = size(x0);
x = abs(x0);

y0 = poissrnd(x);


%%%% generate mask
fractionmissing = 0.75;
q = rand(M,N);
mask = double(q>fractionmissing);

%%%% observation
y = mask.*y0;

figure, imagesc(x0), colormap gray;
figure, imagesc(y0), colormap gray;
figure, imagesc(y), colormap gray;

%%% Estimate signal and mask

% %% add positive bias term to avoid zeros
bias = 1e-3;
%y = y + bias;


ya = Anscombe_forward(y);
figure, imagesc(ya), colormap gray;


y = ya + bias;
g = log(y);
%figure, imagesc(g), colormap gray;


iters = 2500;

% tau1 = 1e-3; % 0.005
% tau2 = ;

%%%% for 0.9
% tau1 = 2.5;
% tau2 = 1;

% tau1 = 6;
% tau2 = 0.5;

%%%% works for 0.9 and 0.7
tau1 = 6;
tau2 = 0.5;


t0 = cputime;
times(1) = 0;

u = zeros(size(g));
v = ones(size(mask))*(-2);


iters = 500;
tau1 = 0.02; % 0.02 for 40% missing
tau2 = 0.1;

u = zeros(size(g));
%v = 01-mask)*(-10);
v = ones(size(g))*(-10);

figure;
x_hat = g;%exp(u);
x_hat_prev = x_hat;

mask_est = exp(v);
mask_est_prev = mask_est;

obj(1) = norm(ya(:)-x_hat(:).*mask_est(:))^2 + tau1*TVnorm(x_hat) + tau2*sum(v(:)~=0);


for iter = 1:iters
    
%     iter
    %u = soft(g-v,tau1);
    %u = chambolledenoise1d(g-v,tau1,100);
    %u = projk(g-v,tau1,5);
    u = tvdenoise(g-v,tau1,5);
    v = hard(g-u,tau2);   
    
    x_hat = exp(u);
    mask_est = exp(v);
    
    times(iter+1) = cputime-t0;
    
    criterion1(iter) = norm(x_hat(:)-x_hat_prev(:))/norm(x_hat(:));
    %criterion2(iter) = norm(mask_est(:)-mask_est_prev(:))/norm(mask_est(:));
    criterion2(iter) = sum(abs(mask(:)-mask_est(:))>=0.2);
    
    obj(iter+1) = norm(ya(:)-x_hat(:).*mask_est(:))^2 + tau1*TVnorm(x_hat) + tau2*sum(v(:)~=0);
    
    if abs(obj(iter+1)-obj(iter))/obj(iter) <= 1e-10
        break;
    end
    
    
    x_hat_prev = x_hat;
    mask_est_prev = mask_est;
    
    x_hat = exp(u);
    imagesc(x_hat), colormap gray;
    drawnow
    
%     tau1 = tau1/0.99;
%     tau2 = tau2/0.99;
end

fprintf('%d iters, %g seconds.\n', iter, times(end))

x_hat = exp(u);
figure, imagesc(exp(v)), colormap gray;
% %a = double(g<=-6);
figure, imagesc(x_hat), colormap gray;

figure, loglog(times, obj)
figure, semilogy(times(2:end), criterion1)
figure, semilogy(times(2:end), criterion2)

x_est = GenAnscombe_inverse_closed_form(x_hat,1e-4);
figure, imagesc(x_est), colormap gray;


% 
% [CA,CD] = dwt(x,'haar');
% z = idwt(CA,CD,'haar');
% 
% %wav = daubcqf(2);
% 
% %[temp1,temp2,lev] = mrdwt(x,wav,L);
% % z = mrdwt_TI2D(x,wav,4);
% % 
% figure, plot(z), hold on, plot(x+2);
