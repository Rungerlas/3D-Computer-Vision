%Read the Teddy stereo pair
teddy_left = double(imread('teddyL.pgm'));
teddy_right = double(imread('teddyR.pgm'));
height = size(teddy_left, 1);
width  = size(teddy_left, 2);
%Define the derivative matrix
matrix = [-1 0 1; -1 0 1; -1 0 1];
%Use matrix to compute Ix and its transpose to compute Iy    
left_Ix = conv2(teddy_left, matrix, 'same');
left_Iy = conv2(teddy_left, matrix', 'same');
right_Ix = conv2(teddy_right, matrix, 'same');
right_Iy = conv2(teddy_right, matrix', 'same');
%Use these first order derivaties to compute I2x, I2y and Ixy at each pixel
left_IxH = left_Ix.*left_Ix;
left_Iy2 = left_Iy.*left_Iy;
right_Ix2 = right_Ix.*right_Ix;
right_Iy2 = right_Iy.*right_Iy;
left_Ixy = left_Ix.*left_Iy;
right_Ixy = right_Ix.*right_Iy;

%Define Gaussian smoothing
gaussian_size = 5;
halfsize = (gaussian_size - 1)/2;
standard_deviation = 1; %set standard deviation of the distribution as 1.0 and kernel size 5¡Á5
[x,y] = meshgrid(-halfsize:halfsize, -halfsize:halfsize);
G_1 = exp(-(x.^2 + y.^2)/(2*standard_deviation^2));
G = G_1./sum(G_1(:)); %average the derivative values in 5X5 windows centered at each pixel
%Apply Gaussian smoothing
left_Ixg = conv2(left_IxH, G, 'same');
left_Iyg = conv2(left_Iy2, G, 'same');
right_Ixg = conv2(right_Ix2, G, 'same');
right_Iyg = conv2(right_Iy2, G, 'same');
left_Ixyg = conv2(left_Ixy, G, 'same');
right_Ixyg = conv2(right_Ixy, G, 'same');

%Compute the Harris operator response function
corner_left = zeros(size(teddy_left));
corner_right = zeros(size(teddy_right));
harris_size = 5;
halfsize2 = (harris_size - 1)/2;
threshold = 98000; %start the value for the threshold from 2000, 
                   %but there will be more than 5000 corners so rise the vaule thershold
                   %until the corners will be 300-500
%The Harris operator  f = det(M)/trace(M)^2
for i = halfsize2+1:height-halfsize2-1
    for j = halfsize2+1:width-halfsize2-1
        %create the matrix M
        left_IxH = sum(sum(left_Ixg(i-halfsize2:i+halfsize2, j-halfsize2:j+halfsize2)));
        left_IyH = sum(sum(left_Iyg(i-halfsize2:i+halfsize2, j-halfsize2:j+halfsize2))); 
        left_IxyH = sum(sum(left_Ixyg(i-halfsize2:i+halfsize2, j-halfsize2:j+halfsize2)));
        right_IxH = sum(sum(right_Ixg(i-halfsize2:i+halfsize2, j-halfsize2:j+halfsize2)));
        right_IyH = sum(sum(right_Iyg(i-halfsize2:i+halfsize2, j-halfsize2:j+halfsize2))); 
        right_IxyH = sum(sum(right_Ixyg(i-halfsize2:i+halfsize2, j-halfsize2:j+halfsize2)));
        %Compare with the threshold, set the corner value and set the rest
        %as 0
        corner_left(i,j) = (left_IxH*left_IyH - left_IxyH*left_IxyH)/(left_IxH+left_IyH);
        %corner_left(i,j) = left_IxyH/(left_IxH^2+left_IyH^2+2*left_IxyH);
        if corner_left(i,j)< threshold
            corner_left(i,j) = 0;
        end            
        corner_right(i,j) = (right_IxH*right_IyH - right_IxyH*right_IxyH)/(right_IxH+right_IyH); 
        %corner_right(i,j) = right_IxyH/(right_IxH^2+right_IyH^2+2*right_IxyH);
        if corner_right(i,j)< threshold
            corner_right(i,j) = 0;
        end  
    end
end
%Apply non-maximum suppression on the responses of the Harris operator in 3X3 windows.
suppression_size = 3;
halfsize3 = (suppression_size - 1)/2;
corner_count1 = 1;
corner_count2 = 1;
for i = halfsize3+1:size(teddy_left,1)-halfsize3-1
  for j = size(teddy_left,2)-halfsize3-1:-1:halfsize3+1
    %judge whether the pixel has the maximum response 
     if corner_left(i,j) ~= 0
            neighborhood = corner_left(i-halfsize3:i+halfsize3, j-halfsize3:j+halfsize3);
            if sum(sum(corner_left(i,j) < neighborhood)) == 0
                corner_left2(corner_count1,:) = [corner_left(i,j) i j]; %save the corners'value and states
                corner_count1 = corner_count1+1;
            else
                corner_left(i,j) = 0; %Make sure that the order processed pixels does
            end                            %not affect the output 
     end        
     if corner_right(i,j) ~= 0
            neighborhood = corner_right(i-halfsize3:i+halfsize3, j-halfsize3:j+halfsize3);
            if sum(sum(corner_right(i,j) < neighborhood)) == 0
                corner_right2(corner_count2,:) = [corner_right(i,j) i j];
                corner_count2 = corner_count2+1;
            else
                corner_right(i,j) = 0;
            end
     end
   end
end
total_left = length(corner_left2);
total_right = length(corner_right2);

%Using SAD to compute the distance
sad_size = 3;
halfsize4 = (sad_size - 1)/2;
corner_count3 = 1;
for i = 1:length(corner_left2)
    for j = 1:length(corner_right2)
       xl = corner_left2(i,2);
       xr = corner_right2(j,2);
       yl = corner_left2(i,3);
       yr = corner_right2(j,3);
       %set 3X3 window for sad
       sad_left = teddy_left(xl-halfsize4:xl+halfsize4, yl-halfsize4:yl+halfsize4);
       sad_right = teddy_right(xr-halfsize4:xr+halfsize4, yr-halfsize4:yr+halfsize4);
       distance(corner_count3,:) = [sum(sum(abs(sad_left-sad_right))) i j];
       corner_count3 = corner_count3+1;
    end
end
%Compute the correct and incorrect
k = 20;
   corner_count5 = 0;
   %get top #% distance
   top_distance = ceil((k/100)*length(distance));
   distance = sortrows(distance);
   distance = distance(1:top_distance,:);
   gt = imread('disp2.pgm');
   gt = double(gt)./4;
   for i = 1:length(distance)
     disparity = distance(i,3);
     disparity_gt = gt(corner_left2(distance(i,2),2), corner_left2(distance(i,2),3));
        if abs(disparity-disparity_gt) <= 1
         corner_count5 = corner_count5 + 1;
        end
   end
disp(['Top' num2str(k) '% most likely correspondences: '])
disp(['Correct number:' num2str(corner_count5)  '   Wrong number: ' num2str(length(distance)-corner_count5)])

