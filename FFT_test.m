clc;
clear;

img = im2double(rgb2gray(imread('lena.jpg')));
img2=img;
img = fft2(img);
img_r = real(img);
img_i = imag(img);
img2_r = hilbert(img_i);
img2_i = -1.*hilbert(img_r);

%  figure(1);
%  subplot(1,2,1);
%  imshow(img_r,[-5,5]);
%  title("real");
%  subplot(1,2,2);
%  imshow(img_i,[-5,5]);
%  title("imag");
%  
%  figure(2);
%  subplot(1,2,1);
%  imshow(log(1+abs(img)),[]);
%  subplot(1,2,2);
%  imshow(angle(img));
 
%  figure(3);
%  subplot(1,3,1);
%  A=abs(fftshift(img));
%  imshow(log(1+abs(fftshift(img))),[]);
%  subplot(1,3,2);
%  Fai=angle(fftshift(img));
%  imshow(angle(fftshift(img)));
% 
% A=ifftshift(A);
% Fai=ifftshift(Fai);
% img_res=A.*exp(Fai.*1i);
% subplot(1,3,3);
% imshow(ifft2(img_res));

figure(4);
subplot(2,2,1);
a_center = f_Cut_for_circleMask(log(1+abs(fftshift(img))),100,100,10,0,1);
imshow(a_center,[]);
a_center = f_Cut_for_circleMask(abs(fftshift(img)),100,100,10,0,1);
a_center = ifftshift(a_center);
Fai=angle(img);
img_res = a_center.*exp(Fai.*1i);
subplot(2,2,2);
imshow(ifft2(img_res));
subplot(2,2,3);
a_out = f_Cut_for_circleMask(log(1+abs(fftshift(img))),100,100,10,0,0);
imshow(a_out,[]);
a_out = f_Cut_for_circleMask(abs(fftshift(img)),100,100,10,0,0);
a_out = ifftshift(a_out);
img_res = a_out.*exp(Fai.*1i);
subplot(2,2,4);
imshow(ifft2(img_res));



function E_out = f_Cut_for_circleMask(E_in,Center_x,Center_y,Cir_Radiu,A,Flag)
% 函数功能：将输入矩阵按照像素值为单位进行圆形区域截取
% E_in：输入矩阵
% Center_x：截取圆形区域的中心坐标x
% Center_y：截取圆形区域的中心坐标x
% Cir_Radiu：截取圆形区域的截取半径像素值
% 注意：这三个参数均是像素值为单位
% A：截取后的其余部分赋值
% Flag:判断截取内部还是外部
[xnums,ynums] = size(E_in);
E_out = E_in;
if Flag
    TAm=1;
    TAn=A;
else
    TAm=A;
    TAn=1;
end
for nx = 1:xnums
    for ny = 1:ynums
        if abs((nx - Center_x) + 1i*(ny - Center_y)) < Cir_Radiu
            E_out(nx, ny) = TAm*E_in(nx, ny);
        else
            E_out(nx, ny) = TAn*E_in(nx, ny);
        end
    end
end
end