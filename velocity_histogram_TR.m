clear all
close all
clc

res =[256 256 256];
spacing = 320/256; %inputs to mimics
frames = 10;
mask_apex = 'Apex_TA.txt';
mask_mid = 'Mid_TA.txt';
mask_base = 'Base_TA.txt';

for frame = 1:frames
    
    frame_name = sprintf('%03d',frame);
    %Read Magnitude Image
    name = sprintf('ph_%03d_mag.dat',frame-1);
    fid = fopen(name,'r');
    raw =fread(fid,'short');
    image = reshape(raw,res);
    fclose(fid);
    
    %Read Magnitude Image
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
    
    %Mimics Read - Apex ------------------------------------------------------
    
    [xM, yM, zM, intensity] = textread(mask_apex,'%f,%f,%f,%f');
    
    xM = (xM/spacing)+1;
    yM = (yM/spacing)+1;
    zM = (zM/spacing)+1;
    
    % disp(['Mimics X Mask Range is ',num2str(min(xM)),' to ',num2str(max(xM))]);
    % disp(['Mimics Y Mask Range is ',num2str(min(yM)),' to ',num2str(max(yM))]);
    % disp(['Mimics Z Mask Range is ',num2str(min(zM)),' to ',num2str(max(zM))]);
    
    v_hist_apex{frame} = zeros(length(xM),1);
    
    for pixel=1:length(xM)
        % xM and yM swapped due to MATLAB indexing protocol
        v_hist_apex{frame}(pixel) = vmag(yM(pixel),xM(pixel),zM(pixel));
    end
     
    %Mimics Read - Mid -------------------------------------------------------
    
    [xM, yM, zM, intensity] = textread(mask_mid,'%f,%f,%f,%f');
    
    xM = (xM/spacing)+1;
    yM = (yM/spacing)+1;
    zM = (zM/spacing)+1;
    
    % disp(['Mimics X Mask Range is ',num2str(min(xM)),' to ',num2str(max(xM))]);
    % disp(['Mimics Y Mask Range is ',num2str(min(yM)),' to ',num2str(max(yM))]);
    % disp(['Mimics Z Mask Range is ',num2str(min(zM)),' to ',num2str(max(zM))]);
    
    v_hist_mid{frame} = zeros(length(xM),1);
    
    for pixel=1:length(xM)
        % xM and yM swapped due to MATLAB indexing protocol
        v_hist_mid{frame}(pixel) = vmag(yM(pixel),xM(pixel),zM(pixel));
    end
    
    %Mimics Read - Base -------------------------------------------------------
    
    [xM, yM, zM, intensity] = textread(mask_base,'%f,%f,%f,%f');
    
    xM = (xM/spacing)+1;
    yM = (yM/spacing)+1;
    zM = (zM/spacing)+1;
    
    % disp(['Mimics X Mask Range is ',num2str(min(xM)),' to ',num2str(max(xM))]);
    % disp(['Mimics Y Mask Range is ',num2str(min(yM)),' to ',num2str(max(yM))]);
    % disp(['Mimics Z Mask Range is ',num2str(min(zM)),' to ',num2str(max(zM))]);
    
    v_hist_base{frame} = zeros(length(xM),1);
    
    for pixel=1:length(xM)
        % xM and yM swapped due to MATLAB indexing protocol
        v_hist_base{frame}(pixel) = vmag(yM(pixel),xM(pixel),zM(pixel));
    end
end

% Plot systole histograms

apex_systole = [v_hist_apex{1}; v_hist_apex{2}; v_hist_apex{3}; v_hist_apex{4}];
mid_systole = [v_hist_mid{1}; v_hist_mid{2}; v_hist_mid{3}; v_hist_mid{4}];
base_systole = [v_hist_base{1}; v_hist_base{2}; v_hist_base{3}; v_hist_base{4}];

figure()
subplot(2,3,1)
histogram(apex_systole,'Normalization','probability','BinWidth',0.005,'EdgeColor','none')
axis([0 1 0 0.1])
xlabel('Velocity [m/s]')
ylabel('Probability density')
title('Apex - Systole')

hold on
avg = mean(apex_systole);
plot(linspace(avg,avg),linspace(0,1),'--')
hold off

disp(['Apex mean: ',num2str(avg)])

subplot(2,3,2)
histogram(mid_systole,'Normalization','probability','BinWidth',0.005,'EdgeColor','none')
axis([0 1 0 0.1])
xlabel('Velocity [m/s]')
ylabel('Probability density')
title('Mid - Systole')

hold on
avg = mean(mid_systole);
plot(linspace(avg,avg),linspace(0,1),'--')
hold off

disp(['Mid mean: ',num2str(avg)])


subplot(2,3,3)
histogram(base_systole,'Normalization','probability','BinWidth',0.005,'EdgeColor','none')
axis([0 1 0 0.1])
xlabel('Velocity [m/s]')
ylabel('Probability density')
title('Base - Systole')

hold on
avg = mean(base_systole);
plot(linspace(avg,avg),linspace(0,1),'--')
hold off

disp(['Base mean: ',num2str(avg)])


% Plot diastole histograms
apex_diastole = [v_hist_apex{6}; v_hist_apex{7}; v_hist_apex{8}; v_hist_apex{9}];
mid_diastole = [v_hist_mid{6}; v_hist_mid{7}; v_hist_mid{8}; v_hist_mid{9}];
base_diastole = [v_hist_base{6}; v_hist_base{7}; v_hist_base{8}; v_hist_base{9}];

subplot(2,3,4)
histogram(apex_diastole,'Normalization','probability','BinWidth',0.005,'EdgeColor','none')
axis([0 1 0 0.1])
xlabel('Velocity [m/s]')
ylabel('Probability density')
title('Apex - Diastole')

hold on
avg = mean(apex_diastole);
plot(linspace(avg,avg),linspace(0,1),'--')
hold off

disp(['Apex mean: ',num2str(avg)])

subplot(2,3,5)
histogram(mid_diastole,'Normalization','probability','BinWidth',0.005,'EdgeColor','none')
axis([0 1 0 0.1])
xlabel('Velocity [m/s]')
ylabel('Probability density')
title('Mid - Diastole')

hold on
avg = mean(mid_diastole);
plot(linspace(avg,avg),linspace(0,1),'--')
hold off

disp(['Mid mean: ',num2str(avg)])


subplot(2,3,6)
histogram(base_diastole,'Normalization','probability','BinWidth',0.005,'EdgeColor','none')
axis([0 1 0 0.1])
xlabel('Velocity [m/s]')
ylabel('Probability density')
title('Base - Diastole')

hold on
avg = mean(base_diastole);
plot(linspace(avg,avg),linspace(0,1),'--')
hold off

disp(['Base mean: ',num2str(avg)])
