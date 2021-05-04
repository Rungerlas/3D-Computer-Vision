clc
clear all
% load variables: BackgroundPointCloudRGB,ForegroundPointCloudRGB,K,crop_region,filter_size)
load data.mat

data3DC = {BackgroundPointCloudRGB,ForegroundPointCloudRGB};
R       = eye(3);
move    = [0 0 -0.25]';
%Get the width and height of the target picture and create a bounding box
%640X400
width = range(ForegroundPointCloudRGB(1,:));
height = range(ForegroundPointCloudRGB(2,:));
z0 = mean(ForegroundPointCloudRGB(3,:)); %to find the Z axis
picture_height = 640;
picture_width = 400;
%Create the file of video for the Dolly Zoom
newvideo = VideoWriter('OutputVideo');
newvideo.FrameRate = 15; %Set the frame rate of 15HZ
open(newvideo); 
%Create the beginnig position of the camera
M0 = K*[R move*0];
sampleposition = [ForegroundPointCloudRGB(1,:);
    ForegroundPointCloudRGB(2,:);
    ForegroundPointCloudRGB(3,:);
    ones(1,size(ForegroundPointCloudRGB(1,:),2))];
new_position = M0*sampleposition;
%To generate the Z and each step moving /Set step as 40
Z = range(new_position(1,:) ./ new_position(3,:));
z1 = z0 * Z / picture_width;
t0 = [0; 0; z1 - z0];
t1 = [0; 0; -3.8];
move = [zeros(2, 75);linspace(0, t1(3)-t0(3), 75)];  %generate new t and new move for the Projection matrix
%Create the capture sequence of image and write into the video file
for step = 1:75
    tic;
    fname       = sprintf('SampleOutput%03d.jpg',step);
    display(sprintf('\nGenerating %s',fname)); 
    K(1,1) = picture_width / width * (z1 + move(3, step)); %Set fx	= (Depthinit + zmove) ¡Á Widthimg¨MWidthobj 
    K(2,2) = picture_height / height * (z1 + move(3, step));%Set fy	= (Depthinit + zmove) ¡Á Heightimg¨MHeightobj
    M = K*[R t0 + move(:, step)];
    im = PointCloud2Image(M, data3DC, crop_region, filter_size);
    imwrite(im, fname);
    writeVideo(newvideo, im);
    toc;
end

close(newvideo);