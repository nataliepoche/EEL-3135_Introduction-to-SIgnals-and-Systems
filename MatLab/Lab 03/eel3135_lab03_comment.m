% USER DEFINED VARIABLES
w = 20;                                                         % Width of the Gaussian function
x = 0:1:79;                                                     % Horiztonal Axis
y = 0:1:79;                                                     % Vertical Axis

% ==> Creates a Gaussian-like 2D matrix based on specified width(w), generated 
% using the exponential function to create a smooth decay from the center 
% with a specified spread <==
z = round(exp(-1/w.^2*(((y.'-50)/1.5).^2+((x-20)).^2)));

% ==> Apply image processing systems to the matrix z.
% The first system (image_system1) samples and scales the image, while the
% second system (image_system2) applies a transformation based on the output
% of the first system <==
[xs,ys,zs] = image_system1(z,3,6);
za         = image_system2(zs,90,-4);

% PLOT RESULT WITH SUBPLOT
figure(1);			
subplot(1,3,1);                                                 % ==> Creates a subplot for the original image <==
imagesc(x, y, z);                                               % ==> Displays the original z image with specified x and y axes <==
axis square; axis xy;	                                        % ==> Sets the aspect ratio to be equal and the y axis to increase up <==
title('Original')                                               % ==> Sets the title of the subplot <==

subplot(1,3,2);                                                 % ==> Creates a subplot for the output of image_system1 <==
imagesc(xs, ys, zs);	                                        % ==> Displays the processed image zs after applying the changes from image_system1 <==
axis square; axis xy;	                                        % ==> Sets the aspect ratio to be equal and the y axis to increase up <==
title('After System 1')

subplot(1,3,3);                                                 % ==> Creates a new subplot for image_system2 <==
imagesc(xs, ys, za);	                                        % ==> Displays the processed image za after applying the changes from image_system2 <==
axis square; axis xy;	                                        % ==> Sets the aspect ratio to be equal and the y axis to increase up <==
title('After System 2')


function [xs, ys, zs] = image_system1(z,Ux,Dy)
%IMAGE_SYSTEM1   ===> This function processes the input image z by sampling and scaling it.
% Inputs:
%   z  - Input image matrix
%   Ux - Upsampling factor in the horizontal direction
%   Dy - Downsampling factor in the vertical direction <===

% ==> Initializes the output matriz zs with zero restricted based on upscaling and downscaling specifications <==
zs = zeros(ceil(size(z,2)/Dy),ceil(Ux*size(z,1))); 

% ==> Create the new vertical axis ys based on Dy <==
ys = 1:ceil(size(z,1)/Dy);
xs = 1:ceil(Ux*size(z,2));                                      % ==> Create the new horizontal axis xs based on the upsampling factor Ux. <==

% ==> Fill the output matrix zs by sampling the input matrix z.
% The input image is downsampled by Dy in the vertical direction and upsampled
% by Ux in the horizontal direction. <==
zs(1:end,1:Ux:end) = z(1:Dy:end,1:end);

end

function [za] = image_system2(z,Sx,Sy)
%IMAGE_SYSTEM2   ===> This function applies a transformation to the input image z.
% Inputs:
%   z  - Input image matrix
%   Sx - Shift in the horizontal direction
%   Sy - Shift in the vertical direction <===

% ====> Initialize the output matrix za with zeros, having the same size as z. <====
za = zeros(size(z,1), size(z,2)); 

for nn = 1:size(z,1)
	for mm = 1:size(z,2)
		% ====> Check if the current pixel (nn, mm) is within the bounds after shifting. <====
		if nn > Sy && nn-Sy < size(z,1) && mm > Sx && mm-Sx < size(z,2)
			% ====> Assign the value from the shifted position in z to za, scaled by a factor of 1/2. <====
			za(nn,mm) = 1/2*z(nn-Sy,mm-Sx);
		end
	end
end


end
