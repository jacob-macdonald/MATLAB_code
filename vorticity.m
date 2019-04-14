clear all
close all

res =[256 256 256];
spacing = 320/256; %inputs to mimics
frames = 15;

% Load in mask data
mask_name = 'RV_grayvalues.txt';
mask = zeros(size(image));
[xM,yM,zM,intensity] = textread(mask_name,'%f,%f,%f,%f');  
xM = (xM/spacing)+1;
yM = abs((yM/spacing)-256);
zM = abs((zM/spacing)-256);

x = 0:spacing:(res(1)*spacing - spacing);
y = 0:spacing:(res(1)*spacing - spacing);
z = 0:spacing:(res(1)*spacing - spacing);

[X,Y,Z] = meshgrid(x,y,z);

vort = zeros(256,256,256,frames);
mask_vort = zeros(length(xM),frames);
vort_norm = zeros(length(xM),frames);
mask_avg = zeros(frames,1);
mask_total = zeros(1,frames);

tic
for frame = 1:frames
   
    %Read complex difference images
    name = sprintf('ph_%03d_vd_1.dat',frame-1);
    fid = fopen(name,'r');
    U =fread(fid,'short');
    U = reshape(U,res);
    fclose(fid);
    
    name = sprintf('ph_%03d_vd_2.dat',frame-1);
    fid = fopen(name,'r');
    V =fread(fid,'short');
    V = reshape(V,res);
    fclose(fid);
    
    name = sprintf('ph_%03d_vd_3.dat',frame-1);
    fid = fopen(name,'r');
    W =fread(fid,'short');
    W = reshape(W,res);
    fclose(fid);
    
    [curlx,curly,curlz,cav] = curl(X,Y,Z,U,V,W);
    
    vort(:,:,:,frame) = sqrt(curlx.^2 + curly.^2 + curlz.^2);
    
    for ii = 1:length(xM)
        mask_vort(ii,frame) = vort(yM(ii),xM(ii),zM(ii),frame);
        vort_norm(ii,frame) = mask_vort(ii,frame) * sqrt(U(yM(ii),xM(ii),zM(ii))^2+V(yM(ii),xM(ii),zM(ii))^2+W(yM(ii),xM(ii),zM(ii))^2);
    end
 
    mask_avg(frame) = mean(mask_vort(:,frame));
    mask_total(frame) = sum(mask_vort(:,frame));
    norm_total(frame) = sum(vort_norm(:,frame));
  
    toc
end

