
% Written by Muhammed Kotan -   (PhD),Sakarya University, Information Systems Engineering Department 
% 
% If you find this codes are usefull please cite following studies
% 
% [1] Kotan, M., Öz, C., & Kahraman, A. (2021). A linearization-based hybrid approach for 3D reconstruction of objects in a single image. 
%        International Journal of Applied Mathematics and Computer Science, 31(3).
% [2] Kotan, M., & Öz, C. (2017). Surface inspection system for industrial components based on shape from shading minimization approach. 
%       Optical Engineering, 56(12), 123105.
%
% In the estimation of illumination directions and implemantation of linearization-based algorithms 
%   we have utilized from below references. 
%    Plese see following references for detailed information.
% 
% 1-    Pentland, A. (1989). Shape information from shading: a theory about human perception. Spatial vision, 4(2-3), 165-182.
% 2-    Ping-Sing, T., & Shah, M. (1994). Shape from shading using linear approximation. Image and Vision computing, 12(8), 487-498.
% 3-    Zheng, Q., & Chellappa, R. (2002). Estimation of illuminant direction, albedo, and shape from shading.
% 4-	Elhabian, S. Y. (2008). Hands on shape from shading. Computer Vision and Image Processing (CVIP) Laboratory, 
%          University of Louisville, Louisville, KY, Technical Report No. SFS08.

function [Z] = Github_HybridLinearSFS()

%Read Image
Img=imread('tent.png');  
if size(Img,3)==3
    Img=rgb2gray(Img);
end 
Img=double(Img);
%Use filter or not (0:on,1:off)
usefilt=0;
filter=3 %default 

%Use Known illumination or not
prompt = 'Estimate illumination direction? Y/N [Y]: ';
str = input(prompt,'s');
if isempty(str)
    str = 'Y';
end

if str=='Y'
   [L] = illumination(Img)
else
    prompt = 'Enter illumination directions(x,y,z. Use ";" between values. Sample : 0.00;0.23;0.043 \n';
   [L] = input(prompt,'s');
end
    
%iterations
it=25; %25 def

%---------------------------------------------- 
% Kotan&Oz
%---------------------------------------------- 

Z=Kotan_Oz(Img,it,usefilt,filter,L);


%---------------------------------------------
%Use timeit to calculate required time ( 0: PM)
%---------------------------------------------
rt=-1

if(rt==0)
    h=@()Kotan_Oz(Img,it,usefilt,filter,L)
    mytime=timeit(h)
end


%--------------------------------------------------
%Plot all figures
%--------------------------------------------------

 
figure(1);
mesh(Z);
title('Kotan & Oz');

figure(2);
imshow(mat2gray(Z));
title('Kotan & Oz');  

end
function [L] =illumination(I)
% See below references.
% [1] Q. Zheng, R. Chellapa, "Estimation of Illuminant Direction, Albedo, and Shape from
%       Shading," IEEE Transactions on Pattern Analysis and Machine Intelligence ,vol. 13, no. 7, pp. 680-
%       702, July, 1991.
% [2] Elhabian, S. Y. (2008). Hands on shape from shading. Computer Vision and Image Processing (CVIP) Laboratory, University of Louisville, 
%       Louisville, KY, Technical Report No. SFS08.
Mu1 = mean(I(:));
Mu2 = mean(mean(I.^2)); 
[Ex,Ey] = gradient(I);
Exy = sqrt(Ex.^2 + Ey.^2);
nEx = Ex ./(Exy + eps); 
nEy = Ey ./(Exy + eps);
avgEx = mean(nEx(:));
avgEy = mean(nEy(:));
gamma = sqrt((6 *(pi^2)* Mu2) - (48 * (Mu1^2)));
%albedo = gamma/pi;
slant = acos((4*Mu1)/gamma);
tilt = atan(avgEy/avgEx);
if tilt < 0
    tilt = tilt + pi;
end
L = [cos(tilt)*sin(slant) sin(tilt)*sin(slant) cos(slant)];
end

function [Z]= Kotan_Oz(I,iteration1,usefilt,filter,L)
    [p,q]=gradient(I);
    Fp = fft2(p);
    Fq = fft2(q);
    [M,N] = size(I);
    ix = L(1)/L(3);
    iy = L(2)/L(3);
     R =(L(3) + Fp .* L(1)+ Fq .* ...
        L(2))./sqrt(1 +Fp.^2 + Fq.^2);
    Z=I-R;
for k = 1 : iteration1
       R =(L(3) + Fp .* L(1)+ Fq .* ...
        L(2))./sqrt(1 + Fp.^2 + Fq.^2);
    R = max(0,R);
    f = I - R;
    df_dZ =(Fp+Fq).*(ix*Fp + iy*Fq + 1)./(sqrt((1 + Fp.^2 + Fq.^2).^3)* ...
        sqrt(1 + ix^2 + iy^2))-(ix+iy)./(sqrt(1 + Fp.^2 + Fq.^2)* ...
        sqrt(1 + ix^2 + iy^2));
    Z = Z - f./(df_dZ + eps); 
    Z_x(2:M,:) = Z(1:M-1,:);
    Z_y(:,2:N) = Z(:,1:N-1);
    Fp = Z - Z_x;
    Fq = Z - Z_y;
end
 
if usefilt==0
Zpm = medfilt2(abs(Z),[filter filter]);
else
    Zpm=abs(Z);
end

Z = flipud(Zpm);
 end
