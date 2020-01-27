function run_FIDIC(varargin)
%RUN_FIDIC   Run fast iterative digital image correlation (FIDIC).
% 
% run_FIDIC(fname_ref,fname_multipage,savename,w0,d0,inc)
% OR
% run_FIDIC(fname_ref,fname_multipage,savename,w0,d0,inc,image_seq)
% 
% This file calls the DIC analysis and runs it for the images in the
% current working directory. Data is saved in mat format in the current
% working directory. 
% 
% REQUIRED INPUTS
% fname_ref:        reference image
%                   Set to empty array [] if there is no separate
%                   reference image
% fname_multipage:  multipage (or single) tif containing current image(s)
% savename:         name to save data
% w0:               target subset size (pix)
% d0:               subset spacing (pix)
% inc:              set to 0 for cumulative; 1 for incremental comparison
% OPTIONAL INPUT
% image_seq:        Vector conataining sequence of images to correlate. If
%                   this is left blank, correlation will occur for all
%                   images.
% 
%                   If no reference image is given, minimum value of im_seq
%                   should be 2, corresponding to the 2nd image in the
%                   multipage stack.
% 
%                   Important note for incremental correlation: The 
%                   sequence indentified by image_seq are the only images
%                   used in the correlation. This means that the
%                   incremental comparison is between the images and in the 
%                   order specified by image_seq.
%
% 
% Main DIC analysis files were written by members of Christian Franck's 
% reserach group. Citation:
%  E. Bar-Kochba, J. Toyjanova, E. Andrews, K.-S. Kim and C. Franck, A fast 
%  iterative digital volume correlation algorithm for large deformations.
%  Experimental Mechanics, 2015, 55, 261–274.
%
% 
% This script written by Jacob Notbohm, University of Wisconsin 2015-2018
% 

% clear;
close all;
clc;

%% --- INPUTS ---

% % --- Inputs called by function ---
% % Uncomment these if you don't want to run as a function
% 
% % Reference image file name. If no reference image (eg, for computing cell
% % velocities), enter []. Note that if [] is entered, num_images must be at
% % least 2.
% % fname_ref = 'slice526_t06.tif'; 
% fname_ref = [];
% % Deformed image filaname (can be a multipage tif)
% fname_multipage = 'pos14.tif'; % Current image (multipage tif)
% % Incremental (inc = 1) or cumulative (inc = 0) correlation
% inc = 1;
% 
% % Name of file of saved data
% savename = 'FIDIC_results_w0=64.mat';
% 
% % Target subset size (pix). Note that smaller subset size is 
% % susceptible to larger noise.
% % This m file assumes equal subset size and spacing, but in general the
% % DIC written by the Franck ground doesn't require equal subset size and
% % spacing.
% w0 = 64;
% 
% % Desired output subset spacing (pix). Typically choose 1/4 of subset size.
% d0 = w0/4;
% 
% % Images to correlate
% image_seq = [];

% Parse function inputs
fname_ref = varargin{1};
fname_multipage = varargin{2};
savename = varargin{3};
w0 = varargin{4};
d0 = varargin{5};
inc = varargin{6};
L = length(varargin);
if L==6
    image_seq = [];
elseif L==7
    image_seq = varargin{7};
else
    error('Wrong number of inputs.');
end

% --- Inputs not called by function ---

% Initial subset size (pix). I recommend 1 or 2 times w0. Other values will
% sometimes work and sometimes give errors.
sSize_initial = 2*w0;


% --- Check inputs ---

if (inc~=0) && (inc~=1)
    error('Variable ''inc'' must be set to either 0 or 1.')
end
if ~isempty(image_seq)
    if isempty(fname_ref) && min(image_seq)==1
        error('If no reference image is given, minimum value of image_seq should be 2.')
    end
end

%% --- LOAD IMAGES ---

% If image_seq isn't an input, then correlate every image
if isempty(image_seq)
    info = imfinfo(fname_multipage);
    num_images = numel(info);
    % Check to see if last image is blank, which occurs for some cases
    imn = imread(fname_multipage,num_images);
    if all(imn(:)==0)
        num_images=num_images-1;
    end
    image_seq = 1:num_images;
else
    % If image_seq is an input, check to make sure it's ascending and warn
    % the user if not
    iseq2 = sort(image_seq,'ascend');
    if image_seq ~= iseq2
        fprintf('Warning: ''image_seq'' is not in ascending order!\n         I hope you know what you''re doing.\n\n');
    end
end

% Load reference image
if ~isempty(fname_ref)
    im1 = imread(fname_ref);
else % If no filename exists, then 1st image of multipage tif is the reference
    im1 = imread(fname_multipage,1);
    % If there is no reference image, make sure smallest value of image_seq 
    % is 2 so that a blank image isn't correlated.
    if image_seq(1)==1
        image_seq = image_seq(2:end);
    end
end

% Load image(s) to correlate
[M, N] = size(im1);
im2 = uint16(zeros(M,N,length(image_seq)));
if length(image_seq)==1
    % If length()==1 then there's only one image to correlate, so only
    % load one image
    im2_k = imread(fname_multipage, image_seq);
    [M2, N2] = size(im2_k); % Need these extra lines to make sure im1 and im2 are same size
    M3=min(M,M2); N3=min(N,N2);
    im2(1:M3,1:N3) = im2_k(1:M3,1:N3);
else
    % Otherwise load all images
    for k=1:length(image_seq)
        im2_k = imread(fname_multipage,image_seq(k));
        [M2, N2] = size(im2_k);
        M3=min(M,M2); N3=min(N,N2);
        im2(1:M3,1:N3,k) = im2_k(1:M3,1:N3);
    end
end

%% --- RUN DIC FOR EACH IMAGE PAIR ---

% Preallocate cell arrays
u_cell = cell(size(im2,3), 1);
v_cell = u_cell;
x_cell = u_cell;
y_cell = u_cell;
c_peak_cell = u_cell;

for k = 1:size(im2,3)
    % k is correlation number, NOT image number. (Only the images used in
    % the correlation were loaded.)
    tic;        % Start timer
    
    % Get reference image
    if inc==1
        if k>1
            imref = im2(:,:,k-1);
        else
            imref = im1;
        end
    elseif inc==0
        imref = im1;
    end
    
    % Get current image
    imcur = im2(:,:,k);
    
    % FIDIC needs images input as a cell array
    I{1} = imref;
    I{2} = imcur;
    % FIDIC needs initial guess for displacements
    U0 = num2cell(zeros(1,2));
    
    % Call FIDIC
    [U, xy, ~, cc, d0] = IDIC(I,[sSize_initial, sSize_initial],U0,d0,w0);
    % INPUTS
    % I: Images
    % sSize_initial: Initial subset size requires 2 components (rows and cols)
    % u0: Initial guess for displacements
    % d0: Desired subset spacing
    % w0: Desired subset size
    % OUTPUTS
    % u is cell array containting displacements
    % xy cell array containing x and y grid points
    % cc is correlation peak
    % d0 is actual final subset spacing (not the same as requested
    % spacing if IDIC didn't converge at a larger spacing)
    
    % Put displacements into cell for saving later. For some reason I get
    % negative values of displacement, so sign is switched below. Partial
    % cause of - sign is that I plot images using standard xy-axes instead
    % of flipping in the y direction as imaging plotting typically does.
    u_cell{k} = -U{1}; % u is a cell array; first entry is x displacement
    v_cell{k} = -U{2};
    c_peak_cell{k} = cc;
    
    tm = toc; % End time of computation
    disp(['Correlation number ',num2str(k),' completed in ',num2str(tm/60),' mins.'])
end

%% --- SAVE DISPLACEMENTS AND PARAMETERS USED TO RUN DIC ---

% Convert from cells to arrays

% Get size of each cell
[MM, NN] = cellfun(@size,u_cell);
% Final output will be the minimum size
M=min(MM(MM>0)); N=min(NN(NN>0));
idx = find(MM==M,1,'first');
% Preallocate
u = zeros(M,N,length(u_cell))*nan;
v = u;
c_peak = u;
% Convert displacements to arrays
for k=1:length(u_cell)
    u_k = u_cell{k};
    v_k = v_cell{k};
    c_k = c_peak_cell{k};
    if numel(u_k)>0
        u(1:M,1:N,k) = u_k(1:M,1:N);
        v(1:M,1:N,k) = v_k(1:M,1:N);
        c_peak(1:M,1:N,k) = c_k(1:M,1:N);
    end
end
% Gridpoints
x = xy{1}(1:M, 1:N);
y = xy{2}(1:M, 1:N);

% The final number of subsets in each direction is larger than the
% size of the image. This is due to zero padding at the edges of the
% image in IDIC.m and DIC.m. Padded data has values of x and y less than 1 
% and greater than the image size.
[M, N] = size(imref);
Ic = x(1,:)>0 & x(1,:)<N;
Ir = y(:,1)>0 & y(:,1)<M;

x = x(Ir,Ic);
y = y(Ir,Ic);
u = u(Ir,Ic,:);
v = v(Ir,Ic,:);
c_peak = c_peak(Ir,Ic);

% Save
save(savename,'u','v','x','y','d0','w0','inc','c_peak');



