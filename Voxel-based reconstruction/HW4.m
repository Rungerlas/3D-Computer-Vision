addpath D:\projects\Matlab\vgg

rawP = [ 776.649963  -298.408539 -32.048386  993.1581875 132.852554  120.885834  -759.210876 1982.174000 0.744869  0.662592  -0.078377 4.629312012;
    431.503540  586.251892  -137.094040 1982.053375 23.799522   1.964373    -657.832764 1725.253500 -0.321776 0.869462  -0.374826 5.538025391;
    -153.607925 722.067139  -127.204468 2182.4950   141.564346  74.195686   -637.070984 1551.185125 -0.769772 0.354474  -0.530847 4.737782227;
    -823.909119 55.557896   -82.577644  2498.20825  -31.429972  42.725830   -777.534546 2083.363250 -0.484634 -0.807611 -0.335998 4.934550781;
    -715.434998 -351.073730 -147.460815 1978.534875 29.429260   -2.156084   -779.121704 2028.892750 0.030776  -0.941587 -0.335361 4.141203125;
    -417.221649 -700.318726 -27.361042  1599.565000 111.925537  -169.101776 -752.020142 1982.983750 0.542421  -0.837170 -0.070180 3.929336426;
    94.934860   -668.213623 -331.895508 769.8633125 -549.403137 -58.174614  -342.555359 1286.971000 0.196630  -0.136065 -0.970991 3.574729736;
    452.159027  -658.943909 -279.703522 883.495000  -262.442566 1.231108    -751.532349 1884.149625 0.776201  0.215114  -0.592653 4.235517090];


P = zeros(3,4,8);

for i=1:8
    for j=1:3
        P(j,1:4,i) = rawP(i,4*(j-1)+1:4*(j-1)+4);
    end
end
%read all silhouettes
silhouette(:,:,1) = double(imread('silh_cam00_00023_0000008550.pbm'));     
silhouette(:,:,2) = double(imread('silh_cam01_00023_0000008550.pbm')); 
silhouette(:,:,3) = double(imread('silh_cam02_00023_0000008550.pbm'));
silhouette(:,:,4) = double(imread('silh_cam03_00023_0000008550.pbm'));
silhouette(:,:,5) = double(imread('silh_cam04_00023_0000008550.pbm'));
silhouette(:,:,6) = double(imread('silh_cam05_00023_0000008550.pbm'));
silhouette(:,:,7) = double(imread('silh_cam06_00023_0000008550.pbm'));
silhouette(:,:,8) = double(imread('silh_cam07_00023_0000008550.pbm'));
%Define a voxel grid with x ranging from -2.5 m to 2.5 m, y from -3 to 3 m and z from 0 to 2.5 m
x_range = [-2.5 2.5];
y_range = [-3 3];
z_range = [0 2.5];
resolution = 0.008; %start with low resolution prior
position_x = x_range(1):resolution:x_range(2);
position_y = y_range(1):resolution:y_range(2);
position_z = z_range(1):resolution:z_range(2);
[x,y,z] = meshgrid(position_x,position_y,position_z);
voxel = [x(:) y(:) z(:)]; %define the number of voxel 
disp(['Total voxel size: ' num2str(length(voxel))]);
%identify surface points
size1 = size(silhouette(:,:,:),1);
size2 = size(silhouette(:,:,:),2);
for i = 1:8
    initial_camera = P(:,:,i);    
    V_pos = initial_camera*[voxel ones(size(voxel(:,1)))]';
    z1 = ceil(V_pos(1,:)./V_pos(3,:));
    z2 = ceil(V_pos(2,:)./V_pos(3,:)); 
    %get the matching voxel from the silhoutte and sign them 
    judge = (z1 > 0) & (z2 > 0) & (z1 <= size2) & (z2 <= size1);
    judge = find(judge);
    %get point range voxel
    z1 = z1(judge);
    z2 = z2(judge); 
    V_pos2 = sub2ind(size(silhouette(:,:,i)),z2,z1);  
    silhouette2 = silhouette(:,:,i);
    np2 = silhouette2(V_pos2) == 1;  
    judge = judge(np2);
    voxel = voxel(judge,:);   %get the point voxel from silhoutte
end
disp(['The final voxel size from silhoutte:' num2str(length(voxel))]);
writeply(voxel,[ ] , 'voxelmodel.ply'); %create ply file to store 
%Try to Coloring the model
%read all camera pictures
image(:,:,:,1) = imread('cam00_00023_0000008550.png');
image(:,:,:,2) = imread('cam01_00023_0000008550.png');
image(:,:,:,3) = imread('cam02_00023_0000008550.png');
image(:,:,:,4) = imread('cam03_00023_0000008550.png');
image(:,:,:,5) = imread('cam04_00023_0000008550.png');
image(:,:,:,6) = imread('cam05_00023_0000008550.png');
image(:,:,:,7) = imread('cam06_00023_0000008550.png');
image(:,:,:,8) = imread('cam07_00023_0000008550.png');
%voxel-base algorithm
for i = 1:size(P,3)
    initial_camera = P(:,:,i); 
    V_pos = initial_camera*[voxel ones(size(voxel(:,1)))]';
    z1 = ceil(V_pos(1,:)./V_pos(3,:));
    z2 = ceil(V_pos(2,:)./V_pos(3,:));    
    judge = sub2ind(size(image(:,:,:,i)),z2,z1,repmat(1,size(z2,1), size(z2,2)));
    image2 = image(:,:,:,i); 
    R(:,i) = image2(judge);    
    judge = sub2ind(size(image2),z2,z1,repmat(2,size(z2,1), size(z2,2)));
    G(:,i) = image2(judge);    
    judge = sub2ind(size(image2),z2,z1,repmat(3,size(z2,1), size(z2,2)));
    B(:,i) = image2(judge);
end
R = mean(R,2);
G = mean(G,2);
B = mean(B,2);
writeply(voxel, [R G B], 'coloring.ply');
 