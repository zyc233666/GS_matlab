clc;
clear all;
close all;
%%
% Parameters
wavelength = 532e-9;  % Wavelength of the propagating wave (in meters)
pixel_size = 10e-6;  % Size of each pixel in the input plane (in meters)
z_prop = 300e-6;  % Propagation distance (in meters)

img = im2double(rgb2gray(imread('lena.jpg')));
input_field = img;
save('input_field.mat','input_field');
% Load the input field
load('input_field.mat', 'input_field');  % Replace 'input_field.mat' with the actual input field file name
slminput_field = double(input_field);  % Convert to double precision if needed

% Size of the input field
[Ny, Nx] = size(slminput_field);
Lx = Nx * pixel_size;
Ly = Ny * pixel_size;

% Frequency grid
dfx = 1 / Lx;
dfy = 1 / Ly;
fx = -Nx/2 : Nx/2-1;
fy = -Ny/2 : Ny/2-1;
fx = fx * dfx;
fy = fy * dfy;
[Fx, Fy] = meshgrid(fx, fy);
kx = 2 * pi * Fx;
ky = 2 * pi * Fy;

% Angular spectrum propagation kernel
prop_kernel = exp(1i * z_prop * sqrt((2 * pi * Fx / Lx).^2 + (2 * pi * Fy / Ly).^2) - 1i * pi * wavelength * z_prop);

% Fourier transform of the input field
input_field_FT = fftshift(fft2(slminput_field));
input_field_cen = f_Cut_for_circleMask(input_field_FT,100,100,10,0,1);
input_field_FT = f_Cut_for_circleMask(input_field_FT,100,100,10,0,0);
% Propagation in the Fourier domain
output_field_FT = input_field_FT .* prop_kernel;
output_field_FT = output_field_FT+input_field_cen;
% Backpropagation to the input plane
output_field = ifft2(ifftshift(output_field_FT));

% Crop the output field to the desired size
Nx_output = round(Lx / pixel_size);
Ny_output = round(Ly / pixel_size);
output_field = output_field(1:Ny_output, 1:Nx_output);
Z_wafer=abs(output_field).^2;
save('Z_wafer.mat','Z_wafer');
% Display the propagated field
figure;
imagesc(abs(output_field).^2);
axis image;
colormap('jet');
colorbar;
title('Propagated Field');
%%
% Parameters
wavelength = 532e-9;  % Wavelength of the propagating wave (in meters)
pixel_size = 10e-6;  % Size of each pixel in the reconstruction plane (in meters)
z_recon = 0.1;  % Reconstruction distance (in meters)
z_prop = 300e-6;  % Propagation distance (in meters)

% Load the input field
load('Z_wafer.mat', 'Z_wafer');  % Replace 'input_field.mat' with the actual input field file name
Z_wafer = double(Z_wafer);  % Convert to double precision if needed

% Size of the input field
[Ny, Nx] = size(Z_wafer);
Lx = Nx * pixel_size;
Ly = Ny * pixel_size;

% Frequency grid
dfx = 1 / Lx;
dfy = 1 / Ly;
fx = -Nx/2 : Nx/2-1;
fy = -Ny/2 : Ny/2-1;
fx = fx * dfx;
fy = fy * dfy;
[Fx, Fy] = meshgrid(fx, fy);
kx = 2 * pi * Fx;
ky = 2 * pi * Fy;

% Angular spectrum propagation kernel
prop_kernel = exp(1i * wavelength * sqrt(1 - (wavelength^2 * (kx.^2 + ky.^2))));
prop_kernel = prop_kernel / (1i * wavelength * z_prop);
prop_kernel = exp(1i * z_prop * sqrt((2 * pi * Fx / Lx).^2 + (2 * pi * Fy / Ly).^2) - 1i * pi * wavelength * z_prop);

% Fourier transform of the input field
input_field_FT = fftshift(fft2(Z_wafer));
input_field_FT = f_Cut_for_circleMask(input_field_FT,100,100,10,0,0);
% Propagation in the Fourier domain
output_field_FT = input_field_FT .* prop_kernel;
output_field_FT = output_field_FT+input_field_cen;
% Backpropagation to the reconstruction plane
output_field = ifft2(ifftshift(output_field_FT));

% Crop the output field to the desired size
Nx_recon = round(Lx / pixel_size);
Ny_recon = round(Ly / pixel_size);
output_field = output_field(1:Ny_recon, 1:Nx_recon);

% Propagate to the reconstruction plane
output_field = output_field * exp(1i * 2 * pi * z_recon / wavelength);

re2=abs(output_field).^2;
load('re1.mat','re1');
err1=sqrt(mean2(abs(re1-img).^2));
err2=sqrt(mean2(abs(re2-img).^2));
disp(err1);
disp(err2);


% Display the reconstructed field
figure;
imagesc(abs(output_field).^2);
% re1=abs(output_field).^2;
% save('re1.mat','re1');
axis image;
colormap('jet');
colorbar;
title('Reconstructed Field');
%%
% Read the entire image
image = imread('1.png');  % Replace 'your_image.jpg' with the actual image file name

% Define the coordinates of the selected area
x_start = 100;  % Starting x-coordinate of the selected area
x_end = 300;  % Ending x-coordinate of the selected area
y_start = 50;  % Starting y-coordinate of the selected area
y_end = 250;  % Ending y-coordinate of the selected area

% Extract the selected area from the image
selected_area = image(y_start:y_end, x_start:x_end, :);

% Display the selected area
figure;
imshow(selected_area);
title('Selected Area');
%%
