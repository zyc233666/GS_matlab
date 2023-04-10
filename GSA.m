clc; clear ; close all;
% based on franuhofer diffraction
RMS_GS = zeros(1000,1);                       %����GS�㷨���������

img = im2double(rgb2gray(imread('lena.jpg')));
img_a = sqrt(img);
disp(class(img));
disp(size(img));
 
% ��slmƽ�濪ʼ
slm_a = 1;
slm_p = rand()*2*pi;
slm_c = img_a.*exp(1i*slm_p);

tic
for i = 1:1000
    % In the image plane, we use the target amplitude
    img_c = fft2(slm_c);
    img_c = img_c./abs(img_c).*img_a; 
    
    % In the slm plane, we use constant amplitude 1 (or gaussian beam
    % amplitude)
    slm_c = ifft2(img_c);
    slm_c = slm_c./abs(slm_c);
    
    % show the orignal target image intensity and its frequency domain
    figure(1);
    subplot(2, 3, 1);
    imshow(img);
    title("Ŀ��ͼ��");
    subplot(2, 3, 2);
    imshow(log(1+abs(fftshift(fft2(img)))), []);
    title("Ŀ��ͼ���Ƶ");
    subplot(2, 3, 3);
    imshow(angle(fftshift(fft2(img))));
    title("Ŀ��ͼ����Ƶ");

    % show the reconstruction image intensity and its frequency domain
    subplot(2, 3, 4);
    reconstruction = power(abs(fft2(slm_c)), 2);
    imshow(reconstruction./max(max(reconstruction)));
    title("�ؽ�ͼ��");
    subplot(2, 3, 5);
    imshow(log(1+abs(fftshift(fft2(reconstruction)))), []);
    title("�ؽ�ͼ���Ƶ");
    subplot(2, 3, 6);
    A=angle(fftshift(fft2(reconstruction)));
    imshow(angle(fftshift(fft2(reconstruction))));
    title("�ؽ�ͼ����Ƶ");
    g_er=abs(img_a)-  fftshift(reconstruction./max(max(abs(reconstruction))));              %�������
    RMS_GS(i)=sqrt(mean2((g_er.^2)));        %������������
%     disp("RMS:"+RMS_GS(i));
    % SSIM > 0.95, then exit the iteration
    if(ssim(img, reconstruction./max(max(reconstruction))) > 0.9)
        N=i;
       disp("Iteration:"+i);
       break; 
    end
end

figure(2)
plot(1:N,RMS_GS(1:N,1));xlabel('ѭ������');ylabel('RMS���(GS)');
% show the orignal target image intensity and its frequency domain
% figure(1);
% subplot(2, 3, 1);
% imshow(img);
% title("Ŀ��ͼ��");
% subplot(2, 3, 2);
% imshow(log(1+abs(fftshift(fft2(img)))), []);
% title("Ŀ��ͼ���Ƶ");
% subplot(2, 3, 3);
% imshow(angle(fftshift(fft2(img))));
% title("Ŀ��ͼ����Ƶ");

% show the reconstruction image intensity and its frequency domain
% subplot(2, 3, 4);
% reconstruction = power(abs(fft2(slm_c)), 2);
% imshow(reconstruction./max(max(reconstruction)));
% title("�ؽ�ͼ��");
% subplot(2, 3, 5);
% imshow(log(1+abs(fftshift(fft2(reconstruction)))), []);
% title("�ؽ�ͼ���Ƶ");
% subplot(2, 3, 6);
% imshow(angle(fftshift(fft2(reconstruction))));
% title("�ؽ�ͼ����Ƶ");
% 
% disp(ssim(img, reconstruction./max(max(reconstruction))));
toc



