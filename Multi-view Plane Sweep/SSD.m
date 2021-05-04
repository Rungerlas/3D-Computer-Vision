function [depth_map] = SSD(camera, x_warp, y_warp, depth, n, Kk, window)
half = floor(window/2);
[rows, cols, color] = size(camera);
depth_map = zeros(rows, cols);
camera = padarray(camera,[half half]);
distance_right = cell(1,length(x_warp));
distance_left = cell(1,length(y_warp));
length1 = length(x_warp);
%Rank Transform
%for i=1:cols
%for j=1:rows
 %define window containing 5X5 pixels around each pixel
 %rankwindow = [max(1,i-half);min(rows,i+half);
 %max(1,j-half);min(cols,j+half)];
 %leftwin=double(x_warp(rankwindow(1):rankwindow(2),rankwindow(3):rankwindow(4)));
 %rightwin=double(y_warp(rankwindow(1):rankwindow(2),rankwindow(3):rankwindow(4)));
 %replace intensity with rank
 %x_warp(i,j) = sum(sum(leftwin<left_img(i,j)));
 %y_warp(i,j) = sum(sum(rightwin<right_img(i,j)));
%end
%end
for disparty = 1:length1
    distance_right{disparty} = (camera-padarray(x_warp{disparty},[half half])).^2;
    distance_left{disparty} = (camera-padarray(y_warp{disparty},[half half])).^2;
end
for j = 1:cols
    for i = 1:rows
        vector = zeros(1, length(x_warp));
        for disparty = 1:length1
            right1 = sum(sum(distance_right{disparty}(i:i+(window-1),j:j+(window-1))));
            left1 = sum(sum(distance_left{disparty}(i:i+(window-1),j:j+(window-1))));
            cost = (right1+left1)/2;
            vector(disparty) = cost;
        end
        [val, index] = min(vector);
        depth_map(i,j) = -depth(index) / ([j i 1.0]*inv(Kk')*n);
    end
end
