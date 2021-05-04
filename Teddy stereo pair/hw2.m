clc; clear all;

%define the range of the disparity
dis_max = 63;
dis_min = 0;
%define the rank window
rank = 5;
window = 3; %3X3 window
%window = 15; %15X15 window
half1 = (rank-1)/2;
half2=(window-1)/2;
%read images
left_img=imread('teddyL.pgm');
right_img=imread('teddyR.pgm');
height=size(left_img,1);
width=size(left_img,2);
%rank transform
ranktleft = zeros(size(left_img));
ranktright =zeros(size(right_img));
 %scan all pixels
for i=1:height
 for j=1:width
   %define window containing 5X5 pixels around each pixel
   rankwindow = [max(1,i-half1);min(height,i+half1);
          max(1,j-half1);min(width,j+half1)]; 
   leftwin=double(left_img(rankwindow(1):rankwindow(2),rankwindow(3):rankwindow(4)));
   rightwin=double(right_img(rankwindow(1):rankwindow(2),rankwindow(3):rankwindow(4)));
   %replace intensity with rank
   ranktleft(i,j) = sum(sum(leftwin<left_img(i,j))); 
   ranktright(i,j) = sum(sum(rightwin<right_img(i,j)));
 end
end
display('Rank transformed finished');

%SAD stereo matching
disparity = zeros(size(left_img));
PKRN = inf(size(left_img)); %set PKRN
 %set cost curve and winner-take
cost=zeros(size(left_img));
curve_entire=inf(size(left_img));
curve_local=inf(size(left_img));
count=1;
%scan all pixels
for i = 1:height
   for j = 1:width
       curve_min=255*5;
       %certain the disparity
        for d = dis_min:dis_max
           cost(i,j)=0;
           %contain pixel window
           if(j+half2-d)>0
           for YL =max(1,i-half2):min(height,i+half2)
             for XL=max(1,j-half2):min(width,j+half2)
                rankleft = ranktleft(YL,XL);
                XR=XL-d; %d=XL-XR
                if XR > 0
                    rankright=ranktright(YL,XL-d);
                end
                cost(i,j)= cost(i,j)+abs(rankleft-rankright); %SAD cost function
             end
           end
           else
               break;
           end
           %compare the cost to find the minima 
           if cost(i,j) <= curve_min
               curve_min = cost(i,j);
               disparity(i,j) = d; %set the pixels' disparity
           end
           %matching cost
           if cost(i,j) < curve_entire(i,j)
               cos = curve_entire(i,j);
               curve_entire(i,j)=cost(i,j);
               curve_local(i,j)=cos;  %the entire curve
           else
               if cost(i,j) <=curve_local(i,j)
                   curve_local(i,j)= cost(i,j); % local minima of the cost curve
               end
           end   
        end
        
        %C PKRN = C2/C1 
        PKRN(i,j)=curve_local(i,j)/curve_entire(i,j);
        if PKRN(i,j) == 1.00
            count=count+1;   %Differences equal to 1 are considered acceptable
        end
    end
end
imshow(disparity);
%imwrite(disparity,'Disparity map(3X3).jpg');
%imwrite(disparity,'Disparity map(15X15).jpg');
 
%error rate
disparity = double(disparity);
gt = imread('disp2.pgm');
gt = double(gt)./4;
errorpoint = sum(sum(abs(gt-disparity)>1))/(size(left_img,1)*size(left_img,2));
display(errorpoint);
imshow(disparity,[dis_min dis_max]);

%computing PKRN error
%generate a disparity map containing the top 50% most con?dent pixels
New_map = sortrows(PKRN);
New_map(1:count) = []; % Pixels without disparity assignments should be ignored in this evaluation
count = 0;
point2 = 0;
for i=1:height
    for j=1:width
        if PKRN(i,j) > New_map(floor(length(New_map)/2))
            
            count = count+1;
        end
    end
end
for i=1:height
    for j=1:width
        if PKRN(i,j) > New_map(floor(length(New_map)/2))
         if abs(gt(i,j)-disparity(i,j))>1
         point2 = point2+1;
         end
        end
    end
end
errorrate = point2/count;
display(errorrate);







