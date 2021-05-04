%Read the old picture and display
picture=imread('basketball-court.ppm');
imshow(picture);
%define new picture size=940X500
new_width=500;
new_height=940;
Output_img = zeros(new_height,new_width,3,'uint8');
%define the target points in the new picture 
%and the points in the old picture
x1 = [1 1];
x2 = [new_width 1];
x3 = [1 new_height];
x4 = [new_width new_height];

xi1 = [247 54];
xi2 = [403 76];
xi3 = [24 195];
xi4 = [281 280];
X = [x1; x2; x3; x4];
Xi = [xi1; xi2; xi3; xi4];
%Create DLT algorithm for homography estimation
%Assemble n 2X9 matrices into a single 2nX9 matrix
Ai=zeros(9,9);
Ai(9,9)=1;
for i = 1:4    
    Ai(2*i-1,1:9)=[0 0 0 -X(i,1) -X(i,2) -1 Xi(i,2)*X(i,1) Xi(i,2)*x(i,2) Xi(i,2)];    
    Ai(2*i,1:9)=[X(i,1) X(i,2) 1 0 0 0 -Xi(i,1)*X(i,1) -Xi(i,1)*X(i,2) -Xi(i,1)];
end   
H0= zeros(9,1); 
H0(9,1)= 1;
H=Ai^-1*H0; 
%Determine h from H
h=reshape(H,[3,3])';
%Create 2D Bilinear interpolation to render the output image
for i=1:new_height
    for j=1:new_width
        point = [(h(1,1)*j+h(1,2)*i+h(1,3))/(h(3,1)*j+h(3,2)*i+h(3,3)) (h(2,1)*j+h(2,2)*i+h(2,3))/(h(3,1)*j+h(3,2)*i+h(3,3))]; %determine the target point 
        xi =floor(point(1)); 
        yj = floor(point(2));
        a = point(1)-xi; %the distance to i
        b = point(2)-yj; %the distance to j
% f(x,y) = (1-a)(1-b)f[i,j]+a(1-b)f[i+1,j]+abf[i+1,j+1]+(1-a)bf[i,j+1]
        Output_img(i,j,:) = picture(yj,xi,:)*(1-a)*(1-b)+picture(yj,xi+1,:)*a*(1-b)+picture(yj+1,xi,:)*b*(1-a)+picture(yj+1,xi+1,:)*a*b;
    end
end
%Display the output image
imshow(Output_img);
imwrite(Output_img,'output1.jpg')