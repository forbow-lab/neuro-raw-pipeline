function [EigenVecs, EigenVals] = kernelEig(kernel, imSize, use_GPU)
% [eigenVecs, eigenVals] = kernelEig(kernel, imSize, [use_GPU=false])
%
% Function computes the ESPIRiT step II -- eigen-value decomposition of a 
% k-space kernel in image space. Kernels should be computed with dat2Kernel
% and then cropped to keep those corresponding to the data space. 
%
% INPUTS:
%           kernel - k-space kernels computed with dat2Kernel (4D)
%           imSize - The size of the image to compute maps for [sx,sy]
%           use_GPU - Flag to use GPU acceleration, default: false
%
% OUTPUTS:
%           EigenVecs - Images representing the Eigenvectors. (sx,sy,Num coils,Num coils)
%           EigenVals - Images representing the EigenValues. (sx,sy,numcoils )
%                       The last are the largest (close to 1)
%           
% 
% See Also:
%               dat2Kernel
% 
%
% (c) Michael Lustig 2010

if(~exist('use_GPU', 'var'));   use_GPU = false;  end
nc = size(kernel,3);
nv = size(kernel,4);
kSize = [size(kernel,1), size(kernel,2)];

% "rotate kernel to order by maximum variance"
kernel = reshape(permute(kernel,[1,2,4,3]),prod([kSize,nv]),nc);

if size(kernel,1) < size(kernel,2)
    [~,~,v] = svd(kernel);
else
    
    [~,~,v] = svd(kernel,'econ');
end

kernel = kernel*v;
kernel = permute(reshape(kernel,[kSize,nv,nc]),[1,2,4,3]);

KERNEL = fft2c(padarray(conj(kernel(end:-1:1,end:-1:1,:,:)), [imSize-kSize,0,0], 'post')) * sqrt(prod(imSize));
KERNEL = permute(KERNEL/sqrt(prod(kSize)), [3,4,1,2]);

if use_GPU
    EigenVecs = zeros(nc, min(nc,nv), prod(imSize));
    EigenVals = zeros(prod(imSize), min(nc,nv));

    for n=1:prod(imSize)
        [EigenVecs(:,:,n), d] = svd(KERNEL(:,:,n),'econ');
        EigenVals (n, :) = real(diag(d));
    end

    EigenVals = reshape(EigenVals(:,end:-1:1), [imSize, min(nc,nv)]);
    
    ph = repmat(exp(-1i*angle(gpuArray(EigenVecs(1,:,:)))),[nc,1,1]);
    EigenVecs = pagefun(@mtimes, v, arrayfun(@(x,y)(x.*y), EigenVecs, ph));
    EigenVecs = gather(reshape(permute(EigenVecs(:,end:-1:1,:), [3, 1, 2]), [imSize, nc, min(nc,nv)]));

else
    EigenVecs = zeros(nc, min(nc,nv), prod(imSize));
    EigenVals = zeros(min(nc,nv), prod(imSize));

    for n=1:prod(imSize)
        [C,D,~] = svd(KERNEL(:,:,n),'econ');

        ph = repmat(exp(-1i*angle(C(1,:))),[size(C,1),1]);
        C = v*(C.*ph);
        D = real(diag(D));
        EigenVals(:,n)   = D(end:-1:1);
        EigenVecs(:,:,n) = C(:,end:-1:1);
    end

    EigenVecs = reshape(permute(EigenVecs, [3, 1, 2]), [imSize, nc, min(nc,nv)]);
    EigenVals = reshape(permute(EigenVals, [2, 1]),    [imSize, min(nc,nv)]);

end
