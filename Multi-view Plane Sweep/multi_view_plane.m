%By reading the ##.png.camera file get the basic data of the picture and
%camera
weight = 3072;
height = 2048;
space = 0.25;
  %The matrices/vectors for rotation, translation and intrinsic parameters for each image
Kk = [2759.48 0 1520.69;0 2764.16 1006.81;0 0 1];
center_rotation = [0.666779 -0.0831384 -0.740603 ;-0.74495 -0.0459057 -0.665539 ;0.021334 0.99548 -0.0925429 ];
center_translation = [-9.46627 -5.58174 0.147736  ]';
left_rotation = [0.582226 -0.0983866 -0.807052  ;-0.813027 -0.0706383 -0.577925 ;-0.000148752 0.992638 -0.121118  ];
left_translation = [-8.31326 -6.3181 0.16107 ]';
right_rotation = [0.795163 -0.050195 -0.604314 ;-0.606377 -0.0736593 -0.791759 ;-0.00477103 0.996019 -0.0890082  ];
right_translation = [-10.8142 -4.53704 0.122293 ]';
%Reset the internal reference of camera and Extend the rotation and translation matrix and reset the reference for the
%homograph
kx = ((weight*space)/2)/(tan(2*atan(Kk(1,3)/Kk(1,1))/2));
ky = ((height*space)/2)/(tan(2*atan(Kk(2,3)/Kk(2,2))/2));
Kk = [kx 0 ((weight*space)/2);0 ky ((height*space)/2);0 0 1];
kk = Kk; %
center_reference = inv(inv([center_rotation center_translation; [0 0 0 1]]));
rotation1 = inv([left_rotation left_translation ; [0 0 0 1]])*center_reference;
rotation2 = inv([right_rotation right_translation ; [0 0 0 1]])*center_reference;
left_rotation = rotation1(1:3,1:3);
left_translation = rotation1(1:3,4);
right_rotation = rotation2(1:3,1:3);
right_translation = rotation2(1:3,4);
%Load the camera image
center_camera = imread('0002.png');
left_camera = imread('0001.png');
right_camera = imread('0003.png');
%Resize the camera for computing
center_camera = imresize(center_camera,space);
left_camera = imresize(left_camera,space);
right_camera = imresize(right_camera,space);
center_camera = im2double(center_camera);
left_camera = im2double(left_camera);
right_camera = im2double(right_camera);
[rows, cols, color] = size(center_camera);
%Set reference for homograph and Cost function
n = right_rotation* [0 0 -1]';
z1 = 5.75;
z2 = 10.75;
plane = 50;
window = 15;
z3 = (z2-z1)/plane;
%Set basic x y d for computing Pixel Dissimilarity
x_warp = cell(1,plane+1);
for x = 1:length(x_warp)
    x_warp{x} = zeros(rows,cols);
end
y_warp = x_warp;
depth = z1:z3:z2;
x1 = repmat(1:cols,rows,1);
y1 = repmat((1:rows)',1,cols);
%Computing homograph and pixels for cost function and matrix P
cost = 0;
for d = z1:z3:z2
    cost = cost + 1;
    Pk1 = kk*(left_rotation-((left_translation*n')/d))*inv(Kk);
    Pk2 = kk*(right_rotation-((right_translation*n')/d))*inv(Kk);
    Pk1 = Pk1/Pk1(3,3);
    Pk2 = Pk2/Pk2(3,3);
    x2 = bsxfun(@plus, bsxfun(@plus, bsxfun(@times, Pk1(1,1), x1), bsxfun(@times, Pk1(1,2), y1)), Pk1(1,3));
    y2 = bsxfun(@plus, bsxfun(@plus, bsxfun(@times, Pk1(2,2), y1), bsxfun(@times, Pk1(2,1), x1)), Pk1(2,3));
    x3 = bsxfun(@plus, bsxfun(@plus, bsxfun(@times, Pk2(1,1), x1), bsxfun(@times, Pk2(1,2), y1)), Pk2(1,3));
    y3 = bsxfun(@plus, bsxfun(@plus, bsxfun(@times, Pk2(2,2), y1), bsxfun(@times, Pk2(2,1), x1)), Pk2(2,3));
    w1  = bsxfun(@plus, bsxfun(@plus, bsxfun(@times, Pk1(3,1), x1), bsxfun(@times, Pk1(3,2), y1)), Pk1(3,3));
    w2  = bsxfun(@plus, bsxfun(@plus, bsxfun(@times, Pk2(3,1), x1), bsxfun(@times, Pk2(3,2), y1)), Pk2(3,3));
    x2 = bsxfun(@rdivide, x2, w1);
    x3 = bsxfun(@rdivide, x3, w2);
    y2 = bsxfun(@rdivide, y2, w1);
    y3 = bsxfun(@rdivide, y3, w2);
    x_warp{cost} = interp2(x1, y1, 255*rgb2gray(left_camera), x2, y2, 'linear', 0);
    y_warp{cost} = interp2(x1, y1, 255*rgb2gray(right_camera), x3, y3, 'linear', 0);      
end
%Apply the Zero-mean Normalized Cross-correlation (NCC) or SAD to generate the
%Depth map
%[depth_map] = SAD(255*rgb2gray(center_camera), x_warp, y_warp, depth, n, Kk, window);
[depth_map] = SSD(255*rgb2gray(center_camera), x_warp, y_warp, depth, n, Kk, window);
%Display the figure of Depth map raw and colored
figure;
title('Depth Map');
imshow(uint8(depth_map*16));
figure;
title('Color Map');
cim = imagesc(depth_map,[z1, z2]);
colormap(jet);
%Figure average pixel error and error map
count = 0;
load data.mat %load the pointcloud data 
background = BackgroundPointCloudRGB(1:3,:);
foreground = ForegroundPointCloudRGB(1:3,:);
pointcloud = [background foreground];
pointcloud(4,:) = 1;
Kk = [[2759.48 0 1520.69;0 2764.16 1006.81;0 0 1] [0 0 0]'];
point = Kk*[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
point_value = point*pointcloud;
point_value(1,:) = point_value(1,:)./point_value(3,:); 
point_value(2,:) = point_value(2,:)./point_value(3,:);
%Generate an error map
generate_map = zeros(point_value(1,end),point_value(2,end));
for i = 1:length(point_value)
    generate_map(round(point_value(2,i)),round(point_value(1,i))) = point_value(3,i);
end
generate_map = im2double(imresize(generate_map,space));
[rows, cols, color] = size(generate_map);
%Computing the average pixel error
error = abs(generate_map - depth_map);
for i = 1:rows
    for j = 1:cols
        if error(i,j) < 1
            count = count+1;
        end
    end
end
average_pixel = (count/(rows*cols));
disp('Average pixel error:');
disp(average_pixel);
figure;
imagesc(error);
colormap jet;

