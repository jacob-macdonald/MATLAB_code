function P=hypr_phantom

P = zeros(256);

Pcl = lower_circle;
Pv = vertical_bars;
Ph = horizontal_bars;
Pcu = upper_circle;

P = Pcl + Pv + Ph + Pcu;

implay(P)

end

function Pcl = lower_circle
    Pcl = zeros(256,256,1,25);
    contrast = ones(1,25) .* 0.2;
    contrast(1) = 1; contrast(2) = 0.9; contrast(3) = 0.8; contrast(4) = 0.7; contrast(5) = 0.6; contrast(6) = 0.5; contrast(7) = 0.4; contrast(8) = 0.3;
    
    radius = 1:20;
    xc = 150;
    yc = 210;
    theta = 0:1:360;

    for i = 1:25
    
    for r = 1:length(radius);
    for t = 1:length(theta);
        x_temp = yc + radius(r)*cosd(theta(t));
        y_temp = xc + radius(r)*sind(theta(t));
        x_l = floor(x_temp);
        x_h = ceil(x_temp);
        y_l = floor(y_temp);
        y_h = ceil(y_temp);
        
        Pcl(x_l,y_l,1,i) = contrast(i);
        Pcl(x_l,y_h,1,i) = contrast(i);
        Pcl(x_h,y_l,1,i) = contrast(i);
        Pcl(x_h,y_h,1,i) = contrast(i);
    end
    end
    end
end

function Pv = vertical_bars

    Pv = zeros(256,256,1,25);
    
    contrast = [.2 .3 .4 .5 .6 .7 .8 .9 1 .9 .8 .7 .6 .5 .4 .3 .2 .2 .2 .2 .2 .2 .2 .2 .2];
    
    for i = 1:25
    Pv(45:236,20:40,1,i) = contrast(i);
    Pv(75:236,55:70,1,i) = contrast(i);
    Pv(95:236,80:90,1,i) = contrast(i);
    Pv(105:236,95:100,1,i) = contrast(i);
    Pv(111:236,103:106,1,i) = contrast(i);
    Pv(114:236,108:109,1,i) =contrast(i);
    end
end

function Ph = horizontal_bars

    Ph = zeros(256,256,1,25);
    contrast = [.2 .2 .2 .2 .2 .2 .2 .2 .2 .3 .4 .5 .6 .7 .8 .9 1 .9 .8 .7 .6 .5 .4 .3 .2];
    
    for i = 1:25
    
    Ph(20:40, 45:236,1,i) = contrast(i);
    Ph(55:70, 75:236,1,i) = contrast(i);
    Ph(80:90, 95:236,1,i) = contrast(i);
    Ph(95:100,105:236,1,i) = contrast(i);
    Ph(103:106,111:236,1,i) = contrast(i);
    Ph(108:109,114:236,1,i) = contrast(i);
    end
end

function Pcu = upper_circle
    Pcu = zeros(256,256,1,25);
    
    contrast = [.2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .3 .4 .5 .6 .7 .8 .9 1];
    
    radius = 1:20;
    xc = 150;
    yc = 210;
    theta = 0:1:360;
    
    for i = 1:25
    for r = 1:length(radius);
    for t = 1:length(theta);
        x_temp = xc + radius(r)*cosd(theta(t));
        y_temp = yc + radius(r)*sind(theta(t));
        x_l = floor(x_temp);
        x_h = ceil(x_temp);
        y_l = floor(y_temp);
        y_h = ceil(y_temp);
        
        Pcu(x_l,y_l,1,i) = contrast(i);
        Pcu(x_l,y_h,1,i) = contrast(i);
        Pcu(x_h,y_l,1,i) = contrast(i);
        Pcu(x_h,y_h,1,i) = contrast(i);
    end
    end
    end
end