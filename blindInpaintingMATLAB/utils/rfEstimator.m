function [atil, btil, y] = rfEstimator(z,opt,randWins);
% *************************************************************************
% Program to estimate brightness, contrast parameters and RF image from log
% compressed ultrasound images
% *************************************************************************
% Jose Seabra // jseabra@isr.ist.utl.pt
% Last update: May 9th 2012 by David Afonso
% *************************************************************************
%
% Instructions:
% [atil, btil, y] = rfEst(z,opt) estimates the RF image y from the input
% log compressed ultrasound image z.
%
% opt = 0, means that all the image is considered for the estimation
% (use this for homogeneous images).
% opt = 1, means that a sliding window is used to compute the parameterers
% which are thereafter averaged (this works for most heterogeneous images).
% opt = 2, let the user select a region where the parameters are estimated
% (this is a better option when the images show dark regions or strong
% specularities)
% atil is the contrast estimator
% btil is the brightness estimator
% y is the RF image
% *************************************************************************
disp('Running RF estimation...');
eps = 1e-3;
switch(opt),
    case(0)
        z(z==0) = [];
        N2 = size(z(:),1);
        r = min(z(:));
        
        atil = sqrt((24/pi^2)*var(z(:)));
        y0 = exp((z-r)./atil) - 1;
        f0 = mean(y0(:).^2/2);
        
        db = 0.01;
        bk = 0:db:r;
        t = exp((r-bk)./atil)-1;
        pk = N2./(atil.*f0).*t.*(t+1).*exp(-N2./(2.*f0).*t.^2);
        
        btil = sum(bk.*pk.*db);
        y = exp((z-btil)/atil) - 1;
        
    case(1)
        [L,M] = size(z);
        r = min(z(:));
        N2 = size(z(:),1);
        
        h = 20;
        va = [];
        vb = [];
        
        for i=1:h/4:L-h,
            for j=1:h/4:M-h,
                
                zc = z(i:i+h,j:j+h);
                
                at = sqrt((24/pi^2)*var(zc(:)));
                y0 = exp((zc-r)./(at+eps)) - 1;
                f0 = mean(y0(:).^2/2);
                
                db = 0.01;
                bk = 0:db:r;
                t = exp((r-bk)./(at+eps))-1;
                pk = N2./(at.*f0+eps).*t.*(t+1).*exp(-N2./(2.*f0+eps).*t.^2);
                bt = sum(bk.*pk.*db);
                va = [va, at];
                vb = [vb, bt];
                
            end
        end
        
        atil = mean(va);
        btil = mean(vb);
        y = exp((z-btil)/atil) - 1;
        
    case(2)
        
        hold on; title('Select homogeneous region (avoid specular and black regions):');
        colormap gray; axis off; axis square;
        [dum,rect]=imcrop;
        I = imcrop(z,rect);
        [L,M] = size(I);
        
        h = 20;
        va = [];
        vb = [];
        
        for i=1:h/4:L-h,
            for j=1:h/4:M-h,
                
                zc = I(i:i+h,j:j+h);
                r = min(zc(:));
                N2 = size(zc(:),1);
                
                at = sqrt((24/pi^2)*var(zc(:)));
                y0 = exp((zc-r)./at) - 1;
                f0 = mean(y0(:).^2/2);
                
                db = 0.01;
                bk = 0:db:r;
                t = exp((r-bk)./at)-1;
                pk = N2./(at.*f0).*t.*(t+1).*exp(-N2./(2.*f0).*t.^2);
                bt = sum(bk.*pk.*db);
                va = [va, at];
                vb = [vb, bt];
                
            end
        end
        
        atil = mean(va);
        btil = mean(vb);
        y = exp((z-btil)/atil) - 1;
        
    case(3)
        
        va = [];
        vb = [];
        for i=1:length(randWins),
            
            zc = randWins(i); zc = zc{1};
            zc(zc==0) = [];
            zc=zc(:);
            
            N2 = size(zc(:),1);
            r = min(zc(:));
            
            at = [];
            bt = [];
            if ~isempty(zc),
                
                aux = sqrt((24/pi^2)*var(zc));
                if aux~=0,
                    at = sqrt((24/pi^2)*var(zc));
                    y0 = exp((zc-r)./at) - 1;
                    f0 = mean(y0.^2/2) + eps;
                    db = 0.01;
                    bk = 0:db:r;
                    t = exp((r-bk)./at)-1;
                    pk = N2./(at.*f0).*t.*(t+1).*exp(-N2./(2.*f0).*t.^2);
                    bt = sum(bk.*pk.*db);
                    
                end
            end
            va = [va, at];
            vb = [vb, bt];
            
            
        end
        
        atil = mean(va);
        btil = mean(vb);
        y = exp((z-btil)/atil) - 1;
    otherwise
        disp('Unknown option in rfEstimatior');
                
end
disp('Done estimating RF.');


end