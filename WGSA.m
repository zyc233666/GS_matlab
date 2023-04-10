%%%%% Copyright 2021.04.05 Beijing %%%%
clc;
clear;
% gray = im2double(rgb2gray(imread('lena.jpg')));
% Amplitude = sqrt(gray);
gray = double(rgb2gray(imread('lena.jpg')));               %�����ʼͼ,��Ϊ��ʼ�����ֲ���� 
Amplitude = imresize(gray,[512,512]);             %��Ӧ�������ֱ��ʣ���ע�͵�
Amplitude = Amplitude./(max(max(Amplitude)));     %��һ��
phase = 2*pi*rand(512,512);                  %���������λ
% phase = rand()*2*pi;
g0_Fie = Amplitude.*exp(1i*phase);           %Fienup�㷨��ʼ������ֲ�
g0_GS = Amplitude.*exp(1i*phase);            %GS�㷨��ʼ������ֲ�
RMS_GS = zeros(500,1);                       %����GS�㷨���������
RMS_Fie = zeros(500,1);                   %����Fienup�㷨���������       



%fienup�㷨,��GS�㷨���뷴��������ger*k;
step_size = 0.8;   %���÷�����������ΧΪ[0,1]��step_size=0ʱΪGS�㷨
for n = 1:500     %�������������� 
%    Fienup�㷨
   G0_Fie = fft2(g0_Fie);            %�渵��Ҷ�任��Ƶ��
   G0_FieNew = 1*G0_Fie./abs(G0_Fie);           %ȡ��λֵ,Ƶ����ȫ1��ֵԼ������λȫϢͼ
   g0_FieNew = ifft2(G0_FieNew);      %������Ҷ�任���ؿ���
   g_er=abs(Amplitude) - abs(g0_FieNew)./max(max(abs(g0_FieNew)));     %������ȷ����������
   RMS_Fie(n)=sqrt(mean2((g_er.^2)));        %������������
   g0_Fie=(abs(Amplitude)+g_er*step_size).*(g0_FieNew./abs(g0_FieNew)); %���뷴������
  
%  GS�㷨
   G0_GS = ifft2(g0_GS);          %�渵��Ҷ�任��Ƶ��
   G0_GSNew = 1*G0_GS./abs(G0_GS);          %ȡ��λֵ,Ƶ����ȫ1��ֵԼ������λȫϢͼ
   g0_GSNew = fft2(G0_GSNew);      %������Ҷ�任���ؿ���
   g_er=abs(Amplitude)-  abs(g0_GSNew)./max(max(abs(g0_GSNew)));              %�������
   RMS_GS(n)=sqrt(mean2((g_er.^2)));        %������������
   g0_GS=abs(Amplitude).*(g0_GSNew./abs(g0_GSNew)); %���뷴������
end

figure(1)
subplot(321);imshow(mat2gray(gray));title('ԭͼ');
subplot(323);imshow(mat2gray(abs(G0_FieNew)));title('��λԭ���ֲ�(Fie)');
subplot(325);imshow(mat2gray(abs(g0_FieNew)));title('ģ���������(Fie)');
subplot(322);imshow(mat2gray(abs(gray)));title('ԭͼ');
subplot(324);imshow(mat2gray(abs(G0_GSNew)));title('��λԭ���ֲ�(GS)');
subplot(326);imshow(mat2gray(abs(g0_GSNew)));title('ģ���������(GS)');

figure(2)
plot(RMS_Fie,'r');
hold on;
plot(RMS_GS,'b');xlabel('ѭ������');ylabel('RMS���(GS)');
lg1=legend('RMS Fie','RMS GS');



