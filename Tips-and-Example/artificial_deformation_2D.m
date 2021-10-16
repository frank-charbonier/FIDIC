function [] = artificial_deformation_2D

clear;
close all;
clc;

% Image to translate
im2 = imread('im02.tif');

% Translations to apply
tx = [2 4 8 0 0 0];
ty = [0 0 0 2 4 8];
% Note that output will be a multipage tif
% image 002 will remain untransformed; it will become slice 1
% slices 002-004 will be translations in x
% slices 005-007 will be translations in y

[M, N] = size(im2);
im_stack = zeros(M, N, length(tx)+1, 'uint16');
im_stack(:,:,1) = im2;

for k=1:length(ty)
    % Apply transformation
    im_stack(:,:,k+1) = tformstack(im2,tx(k),ty(k),0,0);
end

% Save
if exist('translated_images.tif','file')==2
    delete('translated_images.tif');
end
for k=1:size(im_stack,3)
    imwrite(im_stack(:,:,k),'translated_images.tif','writemode','append');
end



%% SUBFUNCTION TO RUN TRANSFORMATION

function transformed_stack = tformstack(im0,tx,ty,ex,ey)
% Apply a computational deformation to an image.
% im0 is the image stack to transform
% tx, ty are translations in column, row directions
% ex, ey are normal strains in y, x directions

% --- Make affine transform ---
T = [(1+ey) 0 ; 0 (1+ex) ; ty tx ];
% Note that y is rows and x is cols (r,c) = (y,x)

tform = maketform('affine',T);

R = makeresampler('cubic', 'fill'); % Specifies type of interpolation and how to handle array boundaries
% These values are more or less the defaults
TDIMS_A = [1 2];
TDIMS_B = [1 2];
TMAP_B = [];
F = [];
TSIZE_B = size(im0);

transformed_stack = tformarray(im0,tform,R,TDIMS_A,TDIMS_B,TSIZE_B,TMAP_B,F);

