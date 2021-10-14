function [u, xy, u_c, cc, dm] = IDIC(varargin)
%ITERATIVE DIGITAL IMAGE CORRELATION
% [u, u_c, cc, dm] = IDVC(I,sSize,u0,dm,wm);
% 
% INPUTS
% -------------------------------------------------------------------------
%   I0: cell containing the undeformed, I0{1}, and deformed, I0{2} images
%   sSize: interrogation window (subset) size
%   u0: pre-estimated displacement field (typically zeros)
%   dm: desired output mesh spacing
%   wm: desired output subset size
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u:   Displacement field vector defined at every meshgrid point with 
%        spacing dm. Format: cell array, each containing a 2D matrix 
%        (components in x,y)
%           u{1} = displacement in x-direction
%           u{2} = displacement in y-direction
%           u{3} = magnitude
%           units: pix
%   xy:  Gridpoints on which displacements u are computed. The gridpoints
%        can vary depending on image size, subset size, and subset spacing.
%        xy is a cell array:
%           xy{1} = grid points giving x positions
%           xy{2} = grid points giving y positions
%           units: pix
%   u_c: same as u but without cropping off some correlations at 
%        zero-padded edges
%   cc:  peak values of the cross-correlation for each interrogation
%   dm:  meshgrid spacing
%
% NOTES
% -------------------------------------------------------------------------
% If used please cite:
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast 
% iterative digital volume correlation algorithm for large deformations. 
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2
% 
% Modified by Jacob Notbohm, University of Wisconsin-Madison, 2016-2018
% 


wb = waitbar(0,'Parsing Inputs','Name','Running IDIC');

% PRESET CONSTANTS
maxIterations = 20; % maximum number of iterations (default 20)
convergenceCrit = [0.25, 0.5, 0.0625]; % convergence criteria
ccThreshold = 1e-4; % bad cross-correlation threshold

[I0, sSize, sSpacing, padSize, DICPadSize, u, dm, wm] = parseInputs(varargin{:});
% sSize and sSpacing are intial subset size and spacing

% START ITERATING
i = 2; % Iteration number is i-1. i==1 corresponds to iteration 0 (ie, the initial subset size and spacing requested).
converged01 = 0; 
SSE = []; I = I0;

t0 = tic;
while ~converged01 && i - 1 < maxIterations
    ti = tic;
    
    set(wb,'name',['Running IDIC (Iteration ',num2str(i-1),')']);
    waitbar(0/7,wb,'Checking Convergence');
    % Check for convergence
     [converged01, SSE(i-1) , sSize(i,:), sSpacing(i,:)] = checkConvergenceSSD_2D(I,SSE,sSize,sSpacing,convergenceCrit,dm,wm);

     % Output subset size and spacing to follow progress
     fprintf('Subset size for iteration %.0f is %.0fx%.0f pix\n',[i-1, sSize(i,:)]);
     fprintf('Subset spacing for iteration %.0f is %.0fx%.0f pix\n',[i-1, sSpacing(i,:)]);
     
    if ~converged01
        [I, m] = parseImages(I,sSize(i,:),sSpacing(i,:));
        
        % run cross-correlation to get an estimate of the displacements
        [du, cc] = DIC(I,sSize(i,:),sSpacing(i,:),DICPadSize,ccThreshold,wm);
  
        % add the displacements from previous iteration to current
        waitbar(3/7,wb,'Adding displacements from previous iteration');
        [u, ~, cc] = addDisplacements_2D(u,du,cc,m,dm);
        
        % filter the  displacements using a predictor filter
        waitbar(4/7,wb,'Filtering Displacements');
        u = filterDisplacements_2D(u,sSize(i,:)/dm);
        
        % remove outliers in displacement field
        waitbar(5/7,wb,'Removing Outliers');
        u = removeOutliers_2D(u);

        % mesh and pad images based on new subset size and spacing
        [I, m] = parseImages(I0,sSize(i,:),sSpacing(i,:));
        
        % map areas based on displacment field
        waitbar(6/7,wb,'Warping Images');
        I = areaMapping_2D(I,m,u);
        
        disp(['Elapsed time (iteration ',num2str(i-1),'): ',num2str(toc(ti))]);
        i = i + 1;
    end
    
end

u_c = u;

% Create 2D set of grid points from mesh grid parameter m
[m{2}, m{1}] = ndgrid(m{1}, m{2});
% Get number of 0s used to pad images
[Itmp, ~] = parseImages(I0,sSize(i-1,:),sSpacing(i-1,:));
[r1, c1] = find(Itmp{1}~=0,1,'first');
xypad = [c1-1, r1-1];
% Function DIC.m pads again by DICPadSize
xypad = xypad + DICPadSize;
% xy grid points are computed by subtracting padding off of mesh grid m
xy{1} = m{1}-xypad(1);
xy{2} = m{2}-xypad(2);
% parseOutputs removes padding
[u, xy, cc] = parseOutputs(u,xy,cc,dm,padSize);

delete(wb);

disp(['Convergence at iteration ',num2str(i)]);
disp(['Total time: ',num2str(toc(t0))]);
end



%% ========================================================================
function varargout = parseImages(varargin)
% pads images and creates meshgrid

I{1} = single(varargin{1}{1});
I{2} = single(varargin{1}{2});
sSize = varargin{2};
sSpacing = varargin{3};

prePad = sSize/2;
postPad = sSize/2;

sizeI = size(I{1});
I{1} = padarray(I{1},prePad,0,'pre');
I{1} = padarray(I{1},postPad,0,'post');

I{2} = padarray(I{2},prePad,0,'pre');
I{2} = padarray(I{2},postPad,0,'post');


idx = cell(1,2);
for i = 1:2, idx{i} = (1:sSpacing(i):(sizeI(i) + 1)) + sSize(i)/2; end

% [m{1},m{2}] = meshgrid(idx{:});

varargout{    1} = I;
varargout{end+1} = idx;

end

%% ========================================================================
function varargout = parseInputs(varargin)
% parses inputs and pads images so that there is an divisable meshgrid number.

I0{1} = single(varargin{1}{1});
I0{2} = single(varargin{1}{2});

% I0{1} = permute(I0{1},[2 1]);
% I0{2} = permute(I0{2},[2 1]);

sSize = varargin{2}; 
sSize = [sSize(2), sSize(1)];

sSpacing = sSize/2;

% Initial guess of displacements
u0 = varargin{3};
% Desired output spacing
dm = varargin{4};
% Desired output subset size
wm = varargin{5};

DICPadSize = sSpacing/2;

sizeI0 = size(I0{1});
sizeI = ceil(sizeI0./sSpacing).*sSpacing;
prePad = ceil((sizeI - sizeI0)/2);
postPad = floor((sizeI - sizeI0)/2);

I{1} = padarray(I0{1},prePad,0,'pre');
I{1} = padarray(I{1},postPad,0,'post');

I{2} = padarray(I0{2},prePad,0,'pre');
I{2} = padarray(I{2},postPad,0,'post');

varargout{    1} = I;
varargout{end+1} = sSize;
varargout{end+1} = sSpacing;
varargout{end+1} = [prePad; postPad];
varargout{end+1} = DICPadSize;
varargout{end+1} = u0;
varargout{end+1} = dm;
varargout{end+1} = wm;
end


function [u, xy, cc] = parseOutputs(u,xy,cc,filterSpacing,padSize)
% parses outputs and unpads the displacment field, gridpoints, and cc.

unpadSize(1,:) = ceil(padSize(1,:)/filterSpacing);
unpadSize(2,:) = floor((padSize(2,:)+1)/filterSpacing); 
% +1 from the extra meshgrid point during the meshing of the DIC algorithm.

for i = 1:2
    % displacements
    u{i} = u{i}(1+unpadSize(1,1):end-unpadSize(2,1),...
        1+unpadSize(1,2):end-unpadSize(2,2));
    % grid points
    xy{i} = xy{i}(1+unpadSize(1,1):end-unpadSize(2,1),...
        1+unpadSize(1,2):end-unpadSize(2,2));
end
u{3} = sqrt(u{1}.^2 + u{2}.^2);

cc = cc(1+unpadSize(1,1):end-unpadSize(2,1),...
    1+unpadSize(1,2):end-unpadSize(2,2));

end

