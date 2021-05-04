function img =PointCloud2Image(M,Sets3DRGB,viewport,filter_size)
 
%% setting up output image
    display('Initializing 2D image...');
    top   = viewport(1);
    left  = viewport(2);
    h     = viewport(3);
    w     = viewport(4);
    bot   = top  + h +1;
    right = left + w +1;
    output_image = zeros(h+1,w+1,3);    
    
    
    for counter = 1:numel(Sets3DRGB)
        display('   Projecting point cloud into image plane...');
        
        % clear drawing area of current layer
        canvas = zeros(bot,right,3);  
        
        % segregate 3D points from color
        dataset          = Sets3DRGB{counter};
        P3D              = dataset(1:3,:);
        color            = dataset(4:6,:)';
        
        % form homogeneous 3D points (4xN)
        X                = [P3D;ones(1,size(P3D,2))];
        
        % apply (3x4) projection matrix
        x                = M*X;
       
        % normalize by 3rd homogeneous coordinate
        x                = x ./ [x(3,:);x(3,:);x(3,:)];
        
        % truncate image coordinates
        x(1:2,:)         = floor(x(1:2,:));

        % determine indices to image points within crop area
        i1     = x(2,:)>top;
        i2     = x(1,:)>left;
        i3     = x(2,:)<bot;
        i4     = x(1,:)<right;
        ix     = i1 & i2 & i3 & i4;
        
        % make reduced copies of image points and corresponding color
        rx     = x(:,ix);
        rcolor = color(ix,:);

        % fill canvas with corresponding color
        for i=1:size(rx,2)    
                canvas(rx(2,i), rx(1,i), : )  = rcolor(i,:);
        end

        %crop canvas to desired output size
        cropped_canvas = canvas(top:top+h,left:left+w,:);
        
  
        
        %filter individual color channel
        display('   Running 2D filters...');
        for i=1:3
           %median filter 
           %img(:,:,i)=medfilt2(img2(:,:,i),filter_size); 
           
           %max filter
           filtered_cropped_canvas(:,:,i)=ordfilt2(cropped_canvas(:,:,i),25,true(5));
           
           %no filter
           %img(:,:,i)=img2(:,:,i);
        end
        
        %get indices of pixel drawn in the current canvas
        drawn_pixels = sum(filtered_cropped_canvas,3);
        idx          = drawn_pixels~=0;
        
        %make a 3-chanel copy of the indices
        idxx(:,:,1)  = idx;
        idxx(:,:,2)  = idx;
        idxx(:,:,3)  = idx;
        
        %erase canvas drawn pixels from the output image
        output_image(idxx)   = 0;
        
        %sum current canvas on top of output image
        output_image         = output_image + filtered_cropped_canvas ;
   
    end
    img = output_image;
    display('Done');
    
end