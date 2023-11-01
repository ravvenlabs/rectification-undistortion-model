                      
%rectification parameters
[leftH, rightH, ~, XBounds, YBounds, ~, intrinsicL, intrinsicR, distCoeffsL, distCoeffsR, displacement, offset] = ...
                    computeRectificationParameters(stereoParams, dimOut, 'valid');
                
                recleftfx=1/intrinsicL(1,1);
                recleftfy=1/intrinsicL(2,2);
                recrightfx=1/intrinsicR(1,1);
                recrightfy=1/intrinsicR(2,2);

           
 
    function [Hl, Hr, Q, xBounds, yBounds, success, Kl, Kr, distCoeffsL, distCoeffsR, displacement, offset] = ...
                computeRectificationParameters(this, imageSize, outputView)
            distCoeffsL=[this.CameraParameters1.RadialDistortion, this.CameraParameters1.TangentialDistortion];
            distCoeffsR=[this.CameraParameters2.RadialDistortion, this.CameraParameters2.TangentialDistortion];
            % Make the two image planes coplanar, by rotating each half way
            [Rl, Rr] = computeHalfRotations(this);
            
            % rotate the transltation vector
            t = Rr * this.TranslationOfCamera2';
            
            % Row align the image planes, by rotating both of them such
            % that the translation vector coinsides with the X-axis.
            RrowAlign = computeRowAlignmentRotation(t);
            
            % combine rotation matrices
            Rrectl = RrowAlign * Rl;
            Rrectr = RrowAlign * Rr;
            
            Kl = this.CameraParameters1.IntrinsicMatrix';
            Kr = this.CameraParameters2.IntrinsicMatrix';
            K_new = computeNewIntrinsics(this);
            
            Hleft  = projective2d((K_new * Rrectl / Kl)');
            Hright = projective2d((K_new * Rrectr / Kr)');
            
            Hl=inv(Hleft.T);
            %H1=fi(H1,1,32,28);
            Hr=inv(Hright.T);
           

            
            % apply row alignment to translation
            t = RrowAlign * t;
            
            [xBounds, yBounds, success] = computeOutputBounds(this, ...
                imageSize, Hleft, Hright, outputView);

             n = xBounds(2) - xBounds(1) + 1;
             m = yBounds(2) - yBounds(1) + 1;
            if n~=imageSize(2) || m~=imageSize(1)
                J=n-imageSize(2);
                if mod(J,2)==1
                    xBounds(2)=xBounds(2)-ceil(J/2);
                    xBounds(1)=xBounds(1)+floor(J/2);
                else
                    xBounds(2)=xBounds(2)-(J/2);
                    xBounds(1)=xBounds(1)+(J/2);
                end
                
                if xBounds(1)<=0
                xBounds(2)=xBounds(2) - xBounds(1) + 1;
                xBounds(1)=1;
                
                end          
                K=m-imageSize(1);
                if mod(K,2)==1
                    yBounds(2)=yBounds(2)-ceil(K/2);
                    yBounds(1)=yBounds(1)+floor(K/2);
                else
                    yBounds(2)=yBounds(2)-(K/2);
                    yBounds(1)=yBounds(1)+(K/2);
                end
                
                if yBounds(1)<=0
                yBounds(2)=yBounds(2) - yBounds(1) + 1;
                yBounds(1)=1;
                end 
            end
          
            % [x, y, disparity, 1] * Q = [X, Y, Z, 1] * w
            cy = K_new(2,3) - yBounds(1);
            cx = K_new(1,3) - xBounds(1);
            f_new = K_new(2,2);
            Q = [1, 0,   0,       -cx;
                0, 1,   0,       -cy;
                0, 0,   0,       f_new;
                0, 0, -1/t(1,1), 0]';
            
             params = this.CameraParameters1;
            
            [dL]=computeMap1(params.IntrinsicMatrix,...
                    params.RadialDistortion, params.TangentialDistortion, ...
                    xBounds, ...
                    yBounds, ...
                    Hleft);
          
          
                
            params = this.CameraParameters2;
               
            [dR]=computeMap1(params.IntrinsicMatrix,...
                    params.RadialDistortion, params.TangentialDistortion, ...
                    xBounds, ...
                    yBounds, ...
                    Hright);
     
               displacement=max(dL, dR);
               z=mod(displacement, 10);
              displacement=displacement-z+20;
              offset=displacement+50;

    end
       
        
     function [displacement]=computeMap1(intrinsicMatrix, radialDist, tangentialDist, XBounds, YBounds, H)
            coder.extrinsic('eval');
            [X, Y] = meshgrid(XBounds(1):XBounds(2),...
                YBounds(1):YBounds(2));
            ptsIn = [X(:) Y(:)]; % remapmex requires singles
%              [X1, Y1] = meshgrid(YBounds(1):YBounds(2),...
%                 XBounds(1):XBounds(2));
%             ptsIn1 = [X1(:) Y1(:)]; % remapmex requires singles
            
            if nargin > 4
                ptsIn = H.transformPointsInverse(ptsIn);
            end
            
            if isempty(coder.target)
%                 ptsOut = distortPoints(ptsIn, ...
%                     intrinsicMatrix', radialDist, tangentialDist);
            
                ptsOut = visionDistortPoints(ptsIn, ...
                    intrinsicMatrix', radialDist, tangentialDist);
            else
                ptsOut = vision.internal.calibration.distortPoints(ptsIn, ...
                    intrinsicMatrix, radialDist, tangentialDist);                
            end
            m = YBounds(2) - YBounds(1) + 1;
            n = XBounds(2) - XBounds(1) + 1;
            map = cast(reshape(ptsOut(:,2),[m n]), 'single') + ...
                    single(0);
            map=ceil(map);  
            dim=size(map);
            diff=zeros(dim(1),1);
            for i=1:m
                j=max(map(i,:));
                k=min(map(i,:));
                diff(i,1)=abs(j-k);
            end
            maxDev=max(diff);
            firstRow=max(map(1,:));


            displacement=maxDev+firstRow;
            eval('clear ptsIn'); % be careful with memory
            
                                             
        end
            %------------------------------------------------------------------
        function [Rl, Rr] = computeHalfRotations(this)
            r = vision.internal.calibration.rodriguesMatrixToVector(this.RotationOfCamera2');
            
            % right half-rotation
            Rr = vision.internal.calibration.rodriguesVectorToMatrix(r / -2);
            
            % left half-rotation
            Rl = Rr';
        end                                
                                    
        function RrowAlign = computeRowAlignmentRotation(t)

                xUnitVector = [1;0;0];
                if dot(xUnitVector, t) < 0
                  xUnitVector = -xUnitVector;
                end

            % find the axis of rotation
            rotationAxis = cross(t,xUnitVector);

            if norm(rotationAxis) == 0 % no rotation
              RrowAlign = eye(3);
            else
           rotationAxis = rotationAxis / norm(rotationAxis);
    
          % find the angle of rotation
          angle = acos(abs(dot(t,xUnitVector))/(norm(t)*norm(xUnitVector)));
    
             rotationAxis = angle * rotationAxis;
    
    % convert the rotation vector into a rotation matrix
        RrowAlign = vision.internal.calibration.rodriguesVectorToMatrix(rotationAxis);
            end
        end
        
        
        
        function K_new = computeNewIntrinsics(this)
            % initialize new camera intrinsics
            Kl = this.CameraParameters1.IntrinsicMatrix';
            Kr = this.CameraParameters2.IntrinsicMatrix';
            
            K_new=Kl;
            
            % find new focal length
            f_new = min([Kr(1,1),Kl(1,1)]);
            
            % set new focal lengths
            K_new(1,1)=f_new; K_new(2,2)=f_new;
            
            % find new y center
            cy_new = (Kr(2,3)+Kl(2,3)) / 2;
            
            % set new y center
            K_new(2,3)= cy_new;
            
            % set the skew to 0
            K_new(1,2) = 0;
        end
        
        


 

function [xBounds, yBounds, success] = computeOutputBounds(this, ...
                        imageSize, Hleft, Hright, outputView)
            
            % find the bounds of the undistorted images
            [xBoundsUndistort1, yBoundsUndistort1] = ...
                computeUndistortBounds(this.CameraParameters1, ...
                imageSize, outputView);
            
            undistortBounds1 = getUndistortCorners(xBoundsUndistort1, yBoundsUndistort1);
            
            [xBoundsUndistort2, yBoundsUndistort2] = ...
                computeUndistortBounds(this.CameraParameters2, ...
                imageSize, outputView);
            undistortBounds2 = getUndistortCorners(xBoundsUndistort2, yBoundsUndistort2);                        
            
            % apply the projective transformation
            outBounds1 =  Hleft.transformPointsForward(undistortBounds1);
            outBounds2 = Hright.transformPointsForward(undistortBounds2);
            
            if strcmp(outputView, 'full')
                 [xBounds, yBounds, success] = computeOutputBoundsFull( ...
                    outBounds1, outBounds2);
            else % valid
                [xBounds, yBounds, success] = computeOutputBoundsValid(...
                    outBounds1, outBounds2);
            end            
        end
                              
   function [xBounds, yBounds] = ...
                computeUndistortBounds(this, imageSize, outputView)          
            if strcmp(outputView, 'same')
                xBounds = [1, imageSize(2)];
                yBounds = [1, imageSize(1)];
            else
                [undistortedMask, xBoundsBig, yBoundsBig] = ...
                    createUndistortedMask(this, imageSize, outputView);
                
                [xBounds, yBounds] = getValidBounds(this, undistortedMask, ...
                        xBoundsBig, yBoundsBig, imageSize);
               
            end
        end
                                     
  function [undistortedMask, xBoundsBig, yBoundsBig] = ...
                createUndistortedMask(this, imageSize, outputView)
                
            % start guessing the undistorted mask with the same size of the
            % original image
            xBounds = [1 imageSize(2)];
            yBounds = [1 imageSize(1)];
            
            [X, Y] = meshgrid(xBounds(1):xBounds(2),yBounds(1):yBounds(2));
            ptsIn = [X(:) Y(:)];
                if isempty(coder.target)
                ptsOut = visionDistortPoints(ptsIn, ...
                    this.IntrinsicMatrix', ...
                    this.RadialDistortion, this.TangentialDistortion);
                else
                ptsOut = vision.internal.calibration.distortPoints(ptsIn, ...
                    this.IntrinsicMatrix, ...
                    this.RadialDistortion, this.TangentialDistortion);                
                end
         
%                 ptsOut = distortPoints(ptsIn, ...
%                     this.IntrinsicMatrix', ...
%                     this.RadialDistortion, this.TangentialDistortion);
                     
            mask = zeros(imageSize, 'uint8');
            
            % each pixel in undistorted image contributes to four pixels in
            % the original image, due to bilinear interpolation
            allPts = [floor(ptsOut); ...
                      floor(ptsOut(:,1)),ceil(ptsOut(:,2)); ...
                      ceil(ptsOut(:,1)),floor(ptsOut(:,2)); ...
                      ceil(ptsOut)];
            insideImage = (allPts(:,1)>=1 & allPts(:,2)>=1 ...
                & allPts(:,1)<=imageSize(2) & allPts(:,2)<=imageSize(1));
            allPts = allPts(insideImage, :);
            indices = sub2ind(imageSize, allPts(:,2), allPts(:,1));
            mask(indices) = 1;            
            numUnmapped = prod(imageSize) - sum(mask(:));
            
            % Grow the output size until all pixels in the original image
            % have been used, or the attempt to grow the output size has
            % failed 5 times when new pixels do not contribute to the
            % mapping.
            if numUnmapped > 0 
                p1 = [xBounds(1), yBounds(1)];
                p2 = [xBounds(2), yBounds(2)];
                numTrials = 0;
                while (numTrials < 5 && numUnmapped > 0)

                    p1 = p1 - 1;
                    p2 = p2 + 1;                    
                    w = p2(1) - p1(1) + 1;
                    h = p2(2) - p1(2) + 1;
                    lastNumUnmapped = numUnmapped;

                    ptsIn = [(p1(1):p1(1)+w-1)', p1(2)*ones(w, 1);...
                             (p1(1):p1(1)+w-1)', p2(2)*ones(w, 1);...
                              p1(1)*ones(h, 1),(p1(2):p1(2)+h-1)';...
                              p2(1)*ones(h, 1),(p1(2):p1(2)+h-1)'];
            
                    if isempty(coder.target)
                            ptsOut = visionDistortPoints(ptsIn, ...
                            this.IntrinsicMatrix', ...
                            this.RadialDistortion, this.TangentialDistortion);
                    else
                            ptsOut = vision.internal.calibration.distortPoints(ptsIn, ...
                            this.IntrinsicMatrix, ...
                            this.RadialDistortion, this.TangentialDistortion);                
                    end                  
%                         ptsOut = distortPoints(ptsIn, ...
%                             this.IntrinsicMatrix', ...
%                             this.RadialDistortion, this.TangentialDistortion);
                 
                    
                    newPts = [floor(ptsOut); ...
                              floor(ptsOut(:,1)),ceil(ptsOut(:,2)); ...
                              ceil(ptsOut(:,1)),floor(ptsOut(:,2)); ...
                              ceil(ptsOut)];
                    insideImage = (newPts(:,1)>=1 & newPts(:,2)>=1 ...
                        & newPts(:,1)<=imageSize(2) & newPts(:,2)<=imageSize(1));
                    newPts = newPts(insideImage, :);
                    indices = sub2ind(imageSize, newPts(:,2), newPts(:,1));
                    mask(indices) = 1;
                    numUnmapped = prod(imageSize) - sum(mask(:));
            
                    if lastNumUnmapped == numUnmapped
                        numTrials = numTrials + 1;
                    else
                        numTrials = 0;
                    end
                                        
                    xBounds = [p1(1), p2(1)];
                    yBounds = [p1(2), p2(2)];
                end
            end
            
            % Compute the mapping with the new output size
            xBoundsBig = xBounds;
            yBoundsBig = yBounds;
            
            mask = ones(imageSize, 'uint8');
            fillValuesMask = cast(0, 'uint8');
            
            myMap = vision.internal.calibration.ImageTransformer;
            myMap.update(mask, this.IntrinsicMatrix, ...
                this.RadialDistortion, this.TangentialDistortion, ...
                outputView, xBoundsBig, yBoundsBig);
            
            undistortedMask = myMap.transformImage(mask, 'nearest', fillValuesMask);
  end  
        
  
  function [xBounds, yBounds] = getValidBounds(this, undistortedMask, ...
                xBoundsBig, yBoundsBig, imageSize)
            
            % Get the boundary
            boundaryPixel = getInitialBoundaryPixel(undistortedMask);
            boundaryPixelsUndistorted = bwtraceboundary(undistortedMask, ...
                boundaryPixel, 'W');
            
            % Convert from R-C to x-y
            boundaryPixelsUndistorted = boundaryPixelsUndistorted(:, [2,1]);
            
            % Convert to the coordinate system of the original image
            boundaryPixelsUndistorted(:, 1) = boundaryPixelsUndistorted(:, 1) + xBoundsBig(1);
            boundaryPixelsUndistorted(:, 2) = boundaryPixelsUndistorted(:, 2) + yBoundsBig(1);
            
            % Apply distortion to turn the boundary back into a rectangle
            boundaryPixelsDistorted = distortPoints(this, boundaryPixelsUndistorted);
            
            % Find the pixels that came from the top, bottom, left, and right edges of
            % the original image.
            tolerance = 7;
            minX = max(1, min(boundaryPixelsDistorted(:, 1)));
            maxX = min(imageSize(2), max(boundaryPixelsDistorted(:, 1)));
            minY = max(1, min(boundaryPixelsDistorted(:, 2)));
            maxY = min(imageSize(1), max(boundaryPixelsDistorted(:, 2)));
            topIdx = abs(boundaryPixelsDistorted(:, 2) - minY) < tolerance;
            botIdx = abs(boundaryPixelsDistorted(:, 2) - maxY) < tolerance;
            leftIdx = abs(boundaryPixelsDistorted(:, 1) - minX) < tolerance;
            rightIdx = abs(boundaryPixelsDistorted(:, 1) - maxX) < tolerance;
                        
            % Find the inscribed rectangle.
            topPixels = boundaryPixelsUndistorted(topIdx, 2);
            botPixels = boundaryPixelsUndistorted(botIdx, 2);
            leftPixels = boundaryPixelsUndistorted(leftIdx, 1);
            rightPixels = boundaryPixelsUndistorted(rightIdx, 1);

            % Check if we can compute the valid bounds at all
            coder.internal.errorIf(isempty(topPixels) || isempty(botPixels) || ...
                isempty(leftPixels) || isempty(rightPixels), ...
                'vision:calibrate:cannotComputeValidBounds');
            
            top = max(topPixels);
            bot = min(botPixels);
            left = max(leftPixels);
            right = min(rightPixels);
            
            % Check if the valid bounds cross
            if isempty(coder.target) && (left > right || top > bot ...
                    || minX > tolerance || maxX < imageSize(2)-tolerance ...
                    || minY > tolerance || maxY < imageSize(1)-tolerance)
                warning(message('vision:calibrate:badValidUndistortBounds'));
            end
            
            xBounds = sort([ceil(left), floor(right)]);
            yBounds = sort([ceil(top), floor(bot)]);
  end
  
  function boundaryPixel = getInitialBoundaryPixel(undistortedMask)

sRow = -1;
sCol = -1;
cx = floor(size(undistortedMask, 2) / 2);
for i = floor(size(undistortedMask, 1)/2):size(undistortedMask, 1)
    if undistortedMask(i, cx) == 0
        sRow = i-1;
        sCol = cx;
        break;
    end
end
if sRow == -1
    sRow = size(undistortedMask, 1);
    sCol = cx;
end
boundaryPixel = [sRow, sCol];
end
        
function [xBounds, yBounds, isValid] = computeOutputBoundsFull(...
    outBounds1, outBounds2)

minXY = min(outBounds1);
maxXY = max(outBounds1);
outBounds1 = [minXY; maxXY];

minXY = min(outBounds2);
maxXY = max(outBounds2);
outBounds2 = [minXY; maxXY];

minXY = round(min([outBounds1(1,:); outBounds2(1,:)]));
maxXY = round(max([outBounds1(2,:); outBounds2(2,:)]));
xBounds = [minXY(1), maxXY(1)];
yBounds = [minXY(2), maxXY(2)];
if minXY(1) >= maxXY(1) || minXY(2) >= maxXY(2)
    isValid = false;
else
    isValid = true;
end
end

function undistortBounds = getUndistortCorners(xBounds, yBounds)
undistortBounds = [xBounds(1), yBounds(1);
    xBounds(2), yBounds(1);
    xBounds(2), yBounds(2);
    xBounds(1), yBounds(2);];
end


                
        function [xBounds, yBounds, isValid] = computeOutputBoundsValid(...
    outBounds1, outBounds2)

% Compute the common rectangular area of the transformed images
outPts = [outBounds1; outBounds2];
xSort   = sort(outPts(:,1));
ySort   = sort(outPts(:,2));
xBounds = zeros(1, 2, 'like', outBounds1);
yBounds = zeros(1, 2, 'like', outBounds2);

outBounds1 = round(outBounds1);
outBounds2 = round(outBounds2);
% Detect if there is a common rectangle area that is large enough
xmin1 = min(outBounds1(:,1));
xmax1 = max(outBounds1(:,1));
xmin2 = min(outBounds2(:,1));
xmax2 = max(outBounds2(:,1));

if (xmin1 >= xmax2) || (xmax1 <= xmin2) % no overlap
    isValid = false;
else
    xBounds(1) = round(xSort(4));
    xBounds(2) = round(xSort(5));
    yBounds(1) = round(ySort(4));
    yBounds(2) = round(ySort(5));
    if xBounds(2)-xBounds(1) < 0.4 * min(xmax1-xmin1, xmax2-xmin2) % not big enough
        isValid = false;
    else
        isValid = true;
    end
end
end

 
 
       
        
        
          
    