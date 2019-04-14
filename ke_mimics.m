clear all
close all

res =[256 256 256];
spacing = 320/256; %inputs to mimics
frames = 15;
rho = 1060;

mask_name = 'LV_grayvalues_Jacob.txt';

%Read Magnitude Image
name = sprintf('MAG.dat');
fid = fopen(name,'r');
raw =fread(fid,'short');
image = reshape(raw,res);
fclose(fid);

for frame = 1:frames

    %Read complex difference images
    name = sprintf('ph_%03d_vd_1.dat',frame-1);
    fid = fopen(name,'r');
    raw =fread(fid,'short');
    vx = reshape(raw,res);
    fclose(fid);
    
    name = sprintf('ph_%03d_vd_2.dat',frame-1);
    fid = fopen(name,'r');
    raw =fread(fid,'short');
    vy = reshape(raw,res);
    fclose(fid);
    
    name = sprintf('ph_%03d_vd_3.dat',frame-1);
    fid = fopen(name,'r');
    raw =fread(fid,'short');
    vz = reshape(raw,res);
    fclose(fid);
    
    vmag = sqrt(vx.*vx + vy.*vy + vz.*vz)/1000;
    
    %Apply mask to velocity
    mask = zeros(size(image));
    [xM,yM,zM,intensity] = textread(mask_name,'%f,%f,%f,%f');
    
    xM = (xM/spacing)+1;
    yM = abs((yM/spacing)-256);
    zM = abs((zM/spacing)-256);
    
    for ii = 1:length(xM)
        mask_velocities(ii) = vmag(yM(ii),xM(ii),zM(ii));
    end
    
    vol = (spacing/1000)^3;
    KE(frame) = 0.5 * rho * vol * sum(mask_velocities.^2);
    
    mask_test = zeros(256,256,256);
    image2 = image;
    for jj = 1:length(xM)
        mask_test(yM(jj),xM(jj),zM(jj)) = 1;
        image2(yM(jj),xM(jj),zM(jj)) = 32000;
    end

end

imtool3D(image2,[])


