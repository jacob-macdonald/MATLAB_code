% This script allows for flow analysis in 2D PC VIPR images. The script
% will import the raw data files, allow the user to designate a
% region-of-interest to integrate blood velocity over, and return plots of the
% time-resolved flow. This script is semi-automated, requiring user input
% only to designate ROIs. 

% Author: Jacob Macdonald
% Date: April 9, 2019

% Input:
%   dir: directory containing reconstructed raw image files. If omitted,
%        will default to current directory.
% Output:
%   flow_series: calculated time-resolved flow over user-defined ROI

function flow_series = flow_analysis_2D(dir)

    % Establish working directory. Default to current directory if working
    % directory is not passed to function.
    curDir = pwd; 
    if nargin<1 || (~ischar(dir) && ~isstring(dir)) || ~exist(dir, 'dir')
        dir = curDir;    
        fprintf('Setting source directory to %s\n', dir);
    end
    dir = char(dir);
    
    % Extract relevant acquisition and reconstruction parameters from image
    % headers.
    fileID = fopen(strcat(dir,'\pcvipr_header.txt'),'r');
    pcvipr_header = textscan(fileID,'%s','Delimiter',' ');
    fclose(fileID);
    
    fileID = fopen(strcat(dir,'\data_header.txt'),'r');
    data_header = textscan(fileID,'%s','Delimiter',' ');
    fclose(fileID);
    
    fov_x = str2num(pcvipr_header{1}{14}); % units: mm
    fov_y = str2num(pcvipr_header{1}{16}); % units: mm
    frames = str2num(pcvipr_header{1}{20});
    matrix_x = str2num(data_header{1}{16});
    matrix_y = str2num(data_header{1}{18});
    
    % Calculating other relevant parameters
    xres = fov_x / matrix_x; % units: mm
    yres = fov_y / matrix_y; % units: mm
    pixel_area = xres * yres; % units: mm^2
    
    % Load in velocity-encoded MRI images for flow quantification
    image = zeros(matrix_y, matrix_x, frames);
    TA_image = zeros(matrix_y, matrix_x);
    
    % For now, we will segment on the time-averaged image. Time-resolved
    % segmentation will be incorporated if needed.
    fileID = fopen(strcat(dir,'\MAG.dat'),'r');
    raw = fread(fileID,'short');
    TA_image(:,:) = reshape(raw,matrix_y,matrix_x);
    fclose(fileID);
    
    % Initialize a round ROI
    dimIM = size(TA_image);
    fig = figure();
    hImage = imshow(TA_image, []);
    enableWL(fig)
    xlabel('Click in vessel of interest | Then press any key')
    title('Define Vessel Mask')
    hp = drawpoint; % Click to define center of ROI
    pause;
    cxy = round(hp.Position);
    delete(hp);
    
    [colImage rowImage] = meshgrid(1:dimIM(1), 1:dimIM(2));
    % Next create the circle in the image.
    radius = 10; % radius of starting image
    BWcir = (rowImage - cxy(2)).^2 + (colImage - cxy(1)).^2 <= radius.^2;

    blocations = bwboundaries(BWcir,'noholes');
    pos = blocations{1};
    pos = fliplr(pos);
    hroi = drawfreehand('Position', pos);
    xlabel('Click & drag to refine ROI | Then press any key')
    pause;
    
    % Extract a binary mask from the user-defined ROI
    tmask = hroi.createMask();
    figure()
    title('Vessel Mask')
    imshow(tmask,[])
    
    % Now we will apply the binary vessel mask to each time frame, extract
    % the velocity value for each pixel in the mask, and multiply by pixel
    % area to calculate blood flow. The net flow over time will be plotted
    % and displayed.
    
    temp_image = zeros(matrix_y, matrix_x);
    flow_curve = zeros(1, frames);
    
    for ii = 0:(frames-1)
        % We only need to acess the velocity file in the z direction
        % (direction 3). Although the reconstruction generates files for
        % the x and y components, these files are empty given the 1D nature
        % of our velocity encoding.
        fileID = fopen(strcat(dir,'\',sprintf('ph_%03d_vd_3.dat',ii)),'r');
        raw = fread(fileID,'short');
        temp_image = reshape(raw,matrix_y,matrix_x);
        fclose(fileID);
        
        velocity_mask = temp_image .* tmask; % velocity units: mm/s
        flow_mask = velocity_mask * pixel_area / 1000; % flow units: mL/s
        flow_curve(ii+1) = -sum(flow_mask,'all');
    end
    
    close all % Clean up extra figures we no longer need
    
    % Finally, let's plot the flow curve.
    figure()
    plot(linspace(1,frames,frames), flow_curve, 'LineWidth',2)
    title('Blood flow over time')
    xlabel('Cardiac Frame')
    ylabel('Flow [mL/s]')
    xlim([1 frames])
    
    flow_series = flow_curve
end