function [u, cc] = DIC(varargin)
% [du, cc] = DIC(I,sSize,sSpacing,ccThreshold) estimates
% displacements between two images through digital image
% correlation.
%
% INPUTS
% -------------------------------------------------------------------------
%   I: cell containing the undeformed, I{1}, and deformed, I{2} 2-D images
%   sSize: interrogation window (subset) size
%   sSpacing: interrogation window (subset) spacing.  Determines window
%             overlap factor
%   DICPadSize: 
%   ccThreshold: threshold value that defines a bad cross-correlation
%   wm: taget subset size
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u: cell containing the displacement field (u{1:2} = {u_x, u_y})
%   cc: peak values of the cross-correlation for each interrogation
%
% NOTES
% -------------------------------------------------------------------------
% all functions are self contained
%
% If used please cite:
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast
% iterative digital volume correlation algorithm for large deformations.
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2
%
% Modified by Jacob Notbohm, University of Wisconsin-Madison, 2018-2021
%

% Parse inputs and create meshgrid
[I,m,mSize,sSize,MTF,M,ccThreshold] = parseInputs(varargin{:});

% Initialize variables
mSize_ = prod(mSize);
u12 = zeros(mSize_,2);
cc = zeros(mSize_,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wb = findall(0,'Tag','TMWWaitbar'); wb = wb(1);

waitbar(1/7,wb,'Estimating Displacements (Time Remaining: )');

for k = 1:mSize_
    
    
    tStart = tic; % begin timer
    %-----------------------------------------------------------------------
    % grab the moving subset from the images
    subst = I{1}(m{1}(k,:),m{2}(k,:));
    B = I{2}(m{1}(k,:),m{2}(k,:));
    
    % multiply by the modular transfer function to alter frequency content
    subst = MTF.*subst;
    B = MTF.*B;
    
    % run cross-correlation
    A = xCorr2(subst,B,sSize); %Custom fft-based correlation - very fast
    
    % find maximum index of the cross-correlaiton
    [cc(k), maxIdx] = max(A(:));
    
    % compute pixel resolution displacements
    [u1, u2] = ind2sub(sSize,maxIdx);
    
    % gather the 3x3 pixel neighborhood around the peak
    try xCorrPeak = reshape(A(u1 + (-1:1), u2 + (-1:1)),9,1);
        % least squares fitting of the peak to calculate sub-pixel displacements
        du12 = lsqPolyFit2(xCorrPeak, M{1}, M{2});
        u12(k,:) = [u1 u2] + du12' - (sSize/2) - 1;
        %-----------------------------------------------------------------------
    catch
        u12(k,:) = nan;
    end
    
    % waitbar calculations (update only every 100 iterations)
    if rem(k,100) == 0
        tRemaining = (toc(tStart)*(mSize_ - k)); % Time remaining for waitbar
        waitbar(1/7*(k/mSize_ + 1),wb,['Estimating Displacements (Time Remaining: ',...
            datestr(datenum(0,0,0,0,0,tRemaining),'MM:SS'),')'])
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reshape displacements and set bad correlations to zero
waitbar(2/7,wb,'Removing Bad Correlations')

cc = reshape(double(cc),mSize);
[cc, ccMask] = removeBadCorrelations(I,cc,ccThreshold);

u{1} = reshape(double(u12(:,2)),mSize).*ccMask;
u{2} = reshape(double(u12(:,1)),mSize).*ccMask;

end

%% ========================================================================
function varargout = parseInputs(varargin)
% Parse inputs and create meshgrid

I{1} = varargin{1}{1};
I{2} = varargin{1}{2};
sSize = varargin{2};
sSpacing = varargin{3};
padSize = varargin{4};
ccThreshold = varargin{5};
wm = varargin{6};

% pad images with zeros so that we don't grab any subset outside of the image
% domain. This would produce an error
I{1} = padarray(I{1},padSize,0,'both');
I{2} = padarray(I{2},padSize,0,'both');
sizeV = size(I{1});

% Initialize Mesh Variables
idx = cell(1,2);
for i = 1:2, idx{i} = (1+padSize(i)) : sSpacing(i) : (sizeV(i)-sSize(i)-padSize(i)+1); end
[m{1},m{2}] = ndgrid(idx{:});


% sSize = [sSize(2) sSize(1)];
mSize = size(m{1});
mSize_ = prod(mSize);

m_ = cell(1,2);
for k = 1:2,
    m_{k} = zeros([mSize_,sSize(k)],'uint16');
    repmat_ = repmat((1:sSize(k))-1,mSize_,1);
    m_{k} = bsxfun(@plus, repmat_,m{k}(:));
end

% Initialize quadratic least squares fitting coefficients
[mx, my] = meshgrid((-1:1),(-1:1));
m = [mx(:), my(:)];

M{1} = zeros(size(m,1),6);
for i = 1:size(m,1)
    x = m(i,1); y = m(i,2);
    M{1}(i,:) = [1,x,y,x^2,x*y,y^2];
end

M{2} = M{1}'*M{1};

% Generate Moduluar transfer function (see eq. 3). The 2nd input is the
% equation number used to produce the MTF from the paper cited in the
% subfunction generateMTF. Options are 4, 5, or 6. See comments of 
% generateMTF for more info.
MTF = generateMTF(sSize, wm, 6);

%% Parse outputs

varargout{    1} = I;
varargout{end+1} = m_;
varargout{end+1} = mSize;
varargout{end+1} = sSize;
varargout{end+1} = MTF;
varargout{end+1} = M;
varargout{end+1} = ccThreshold;

end

%% ========================================================================
function A = xCorr2(A,B,sSize)
% performs fft based cross correlation of A and B (see equation 2)

A = fftn(A,sSize);
B = fftn(B,sSize);
B = conj(B);
A = A.*B;
A = ifftn(A);
A = real(A);
A = fftshift(A);
end

%% ========================================================================
function    duvw = lsqPolyFit2(b, M, trMM)
% LeastSqPoly performs a 2D polynomial fit in the least squares sense
% Solves M*x = b,
% trMM = transpose(M)*M
% trMb = tranpose(M)*b
%
% If you need to generate the coefficients then uncomment the following
% [mx, my, mz] = meshgrid(-1:1,-1:1,-1:1);
% m = [mx(:), my(:), mz(:)];
%
% for i = 1:size(m,1)
%    x = m(i,1); y = m(i,2); z = m(i,3);
%    M1(i,:) = [1,x,y,z,x^2,x*y,x*z,y^2,y*z,z^2];  %3D case
%               1 2 3 4  5   6   7   8   9   10

%    M1(i,:) = [1,x,y,x^2,x*y,y^2];                %2D Case
%               1 2 3  5   6   8
%               1 2 3  4   5   6
% end
%
% trMM1 = M'*M;

% b = log(b);
trMb = sum(bsxfun(@times, M, b));

x = trMM\trMb'; %solve for unknown coefficients

A = [x(5), 2*x(4);
    2*x(6),  x(5)];

duvw = (A\(-x([2 3])));

end

%% ========================================================================
function MTF = generateMTF(sSize, wm, eq_no)
% MTF functions taken from
% J. Nogueira, A Lecuona, P. A. Rodriguez, J. A. Alfaro, and A. Acosta.
% Limits on the resolution of correlation PIV iterative methods. Practical
% implementation and design of weighting functions. Exp. Fluids,
% 39(2):314{321, July 2005. doi: 10.1007/s00348-005-1017-1

% Inputs
% - sSize is actual subset size
% - wm is target subset size
% - The option eq_no specifies which equation for the MTF is used. Options
% are taken from different equations in Nogueira et al. Suitable options
% are 4, 5, or 6 corresponding to Eqs. 4, 5, and 6 of the paper. Note that
% the manuscript suggests Eq. 6 as the ``best'' choice.

% Choose condition to generate the MTF. If the condition below isn't
% satisfied, the MTF used is just a constant array and hence has no effect
% if prod(single(sSize ~= 32)) % original by franck group
if prod(single(sSize ~= wm)) % JN suggestion: apply MTF when the target subset is reached; before target is reached, the correlation is probably doing large scale drift correction, in which case using the MTF could cause a correlation to a local maximum rather than the globl max (only works for square subsets)
    
    nu = ones(sSize(1),sSize(2));
    
else
    
    % MTF appears to be written assuming square subsets
    sSize = sSize(1);
    
    switch eq_no
        
        %% equation 4
        case 4
            
            x = cell(1,3);
            [x{1}, x{2}] = meshgrid(1:sSize,1:sSize);
            
            nu = 1;
            for i = 1:2
                x{i} = x{i} - sSize/2 - 0.5;
                x{i} = abs(x{i}/sSize);
                nu = nu.*(3*(4*x{i}.^2-4*x{i}+1));
            end
            
        %% equation 5
        case 5
            
            [x{1}, x{2}] = meshgrid(1:sSize,1:sSize);
            
            for i = 1:2
                x{i} = x{i} - sSize/2 - 0.5; 
            end
            
            r = abs(sqrt(x{1}.^2 + x{2}.^2)/sSize);
            nu  = zeros(size(r));
            nu(r < 0.5) = 24/pi*(4*r(r < 0.5).^2-4*r(r < 0.5)+1);
            
        %% equation 6
        case 6
            [x{1}, x{2}] = meshgrid(1:sSize,1:sSize);
            
            nu = 1;
            for i = 1:2
                x{i} = x{i} - sSize/2 - 0.5;
                x{i} = (x{i}/sSize);
                
                nu = nu.*(12*abs(x{i}).^2 - 12*abs(x{i}) + 3 + ...
                    0.15*cos(4*pi*x{i}) + 0.20*cos(6*pi*x{i}) + ...
                    0.10*cos(8*pi*x{i}) + 0.05*cos(10*pi*x{i}));
            end
            nu(nu < 0) = 0;
            
    end
    
        
    % all equations coded above are actually for the square of nu, so take sqrt
    nu = sqrt(nu);
    
%     % normalize MTF
%     nu = nu/sum(nu(:));
    
    
%     nu = cellfun(@(x) x/sum(x(:)), nu, 'UniformOutput',0);
%     nu = cellfun(@sqrt, nu, 'UniformOutput',0);
    
end

% Output MTF
MTF = nu;

end

%% ========================================================================
function [cc, ccMask] = removeBadCorrelations(I,cc,ccThreshold)
% removes bad correlations.  You can insert your own method here.
minOS = 1;
for i = 1:2
    zeroIdx = I{i} == 0;
    threshold = mean2(I{i}(~zeroIdx));
    I_ = I{i}.*(I{i} < threshold);
    minOS = minOS*sum(I_(:))/sum(~zeroIdx(:));
end
cc = cc - minOS;
cc = cc/(max(I{1}(:))*max(I{2}(:)));
ccMask = double(cc >= ccThreshold);

if prod(ccMask(:)) == 0 % IF statement added 4/12/18 BY JN. This avoids an error when prod(ccMask(:))==1
    CC = bwconncomp(~ccMask);
    [~,idx] = max(cellfun(@numel,CC.PixelIdxList));
    ccMask(CC.PixelIdxList{idx}) = inf;
    ccMask(cc == 0) = nan;
    ccMask(~isfinite(ccMask)) = 0;
end

cc = cc.*ccMask;

end