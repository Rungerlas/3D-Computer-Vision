function [depth_map] = SAD(camera_reference, x_warp, y_warp, depth, n, Kk, window)
half_window = floor(window/2);
[rows, cols, color] = size(camera_reference);
depth_map = zeros(rows, cols);
camera_reference = padarray(camera_reference,[half_window half_window]);
right = cell(1,length(x_warp));
left = cell(1,length(y_warp));
length1 = length(x_warp);
%Apply Rank Transform
%for i=1:height
 %for j=1:width
   %define window containing 5X5 pixels around each pixel
  % rankwindow = [max(1,i-half1);min(height,i+half1);
   %       max(1,j-half1);min(width,j+half1)]; 
   %leftwin=double(left_img(rankwindow(1):rankwindow(2),rankwindow(3):rankwindow(4)));
   %rightwin=double(right_img(rankwindow(1):rankwindow(2),rankwindow(3):rankwindow(4)));
   %replace intensity with rank
   %ranktleft(i,j) = sum(sum(leftwin<left_img(i,j))); 
   %ranktright(i,j) = sum(sum(rightwin<right_img(i,j)));
 %end
%end
for d = 1:length1
    right{d} = abs(camera_reference-padarray(x_warp{d},[half_window half_window]));
    left{d} = abs(camera_reference-padarray(y_warp{d},[half_window half_window]));
end
for j = 1:cols
    for i = 1:rows
        vector = zeros(1, length(x_warp));
        for d = 1:length(x_warp)
            right_cost = sum(sum(right{d}(i:i+(window-1),j:j+(window-1))));
            left_cost = sum(sum(left{d}(i:i+(window-1),j:j+(window-1))));
            cost = (right_cost+left_cost)/2;
            vector(d) = cost;
        end
        [val, index] = min(vector);
        depth_map(i,j) = -depth(index) / ([j i 1.0]*inv(Kk')*n);
    end
end
