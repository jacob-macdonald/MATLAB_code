function [HYPR,FBP] = HYPR_recon(npr,inter)
%   DESCRIPTION: Implementation of the HYPR reconstruction algorithm 
%       described by Wieben et al.
%   INPUTS:
%       npr = number of projections per interleave
%       inter = number of interleaves (time frames)
%   OUTPUTS:
%       HYPR = HYPR reconstructed time-frames
%       FBP = corresponding filtered backprojection time-freams
%   AUTHOR: Jacob Macdonald; November 21, 2014

nproj = npr/inter;

% Create circle phantom (512x512x15)
% P = circ_phantom;
% P = reshape(P,[512 512 inter]);

% Create Shepp-Logan phantom with noise and designate SNR ROI's
P = dynamic_Shepp;
P_SNR = P(76:437,76:437,:);
img = imshow(P_SNR(:,:,1));
title('Designate ROI of signal region')
BW = roipoly;
signal_index = find(BW);
close;

img = imshow(P_SNR(:,:,1));
title('Designate ROI of background noise')
BW2 = roipoly;
noise_index = find(BW2);
close;

SNR_Shepp = linspace(0.3,1,29)/0.05;

% Create interleaved, evenly spaced radial profiles
[kx,ky] = radial_trajectory_interleaved(npr,512,1,512,inter);
kx = fix(kx)+256;
ky = fix(ky)+256;

% Obtain corresponding angles for each projection
theta = zeros(nproj,inter);
angle = linspace(0,180-180/npr,npr);
for ii = 1:inter
    for jj = 1:nproj
        theta(jj,ii) = angle(ii+(inter*jj)-inter);
    end
end

% Obtain projection values in k-space
P_k = fftshift(fftshift(fft2(fftshift(fftshift(P,1),2)),1),2);
P_k_proj = zeros(512,nproj,inter);

for ii = 1:512
    for jj = 1:nproj
        for kk = 1:inter
            P_k_proj(ii,jj,kk) = P_k(ky(ii,jj,kk),kx(ii,jj,kk),kk);
        end
    end
end

% Reshape into composite projections
P_k_proj_comp = zeros(512,npr);

for ii = 1:inter
    P_k_proj_comp(:,(nproj*ii-(nproj-1)):(nproj*ii)) = P_k_proj(:,:,ii);
end

theta_comp = reshape(theta,[1 npr]);

% Take inverse 1D-FFT to obtain projections in image space
P_proj_comp = fftshift(ifft(fftshift(P_k_proj_comp,1)),1);

for ii = 2:inter
    P_comp(:,:,ii-1) = iradon(P_proj_comp(:,1:(nproj*ii)),theta_comp(1,1:(nproj*ii)));
end

% Obtain HYPR time frames
HYPR = zeros(362,362,inter-1);
FBP = zeros(362,362,inter-1);

writerObj = VideoWriter('shepp_hypr_87.avi');
writerObj.FrameRate = 10;
open(writerObj);
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

for ii = 2:inter
    sino_proj = fftshift(ifft(fftshift(P_k_proj(:,:,ii),1)),1);
    phi = theta(:,ii)';
    FBP(:,:,ii-1) = iradon(sino_proj,phi);
    sino_comp = radon(P_comp(:,:,ii-1),phi);
    sino_comp = sino_comp(3:514,:);
    sino_norm = mat2gray(abs(sino_proj)) ./ mat2gray(sino_comp);
    w(:,:,ii-1) = iradon(sino_norm,phi,'Linear','None');
    HYPR_temp = w(:,:,ii-1) .* P_comp(:,:,ii-1);
    TF = isinf(HYPR_temp);
    index = find(TF);
    HYPR_temp(index) = 0;
    HYPR_temp = flipdim(HYPR_temp,1);
    
    SNR_HYPR(ii) = mean(HYPR_temp(signal_index))/std(HYPR_temp(noise_index));
    
    HYPR(:,:,ii-1) = HYPR_temp;
    imshow(HYPR(:,:,ii-1));
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end

close(writerObj);

% Obtain filtered backprojection images

writerObj = VideoWriter('shepp_fbp_87.avi');
writerObj.FrameRate = 10;
open(writerObj);
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

for ii = 2:inter
    sino_proj = fftshift(ifft(fftshift(P_k_proj(:,:,ii),1)),1);
    phi = theta(:,ii)';
    FBP_temp = iradon(sino_proj,phi);
    FBP_temp = flipdim(FBP_temp,1);
    
    SNR_FBP(ii) = mean(FBP_temp(signal_index))/std(FBP_temp(noise_index));
    
    FBP(:,:,ii-1) = FBP_temp;
    imshow(FBP(:,:,ii-1));
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end

close(writerObj);

frame = 1:29;

plot(frame, SNR_HYPR,'*',frame,SNR_FBP,'*');
xlabel('Time frame')
ylabel('SNR estimate')
title('SNR Comparison for HYPR and Filtered Backprojection')
legend('HYPR','FBP','Location','northwest')

end

