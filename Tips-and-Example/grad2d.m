function [dudx, dudy] = grad2d(varargin)
% 
% [dudx, dudy] = grad2d(u)
% 
% OR
% 
% [dudx, dudy] = grad2d(u, n)
% 
% Numerical 2D displacement gradient. Computed by fitting a n x n window of
% points to a line and taking the slope. For computational efficiency, this
% is done with a convolution.
% 
% INPUTS    u: displacement component
%           n: (optional) Size of window. n must be odd. If not specified,
%              the default is a 3x3 window.
% 
% This code assumes unit spacing. If spacing is not 1, divide outputs by 
% spacing in x and y directions, respectively.
% 
% Written by Christian Franck and Soonsung Hong, California Institute of
% Technology 2007.
% 
% Updated by Jacob Notbohm, University of Wisconsin-Madison, 2016-2017
% 

% Inputs
if length(varargin)==1
    u = varargin{1};
    n = 3;
elseif length(varargin)==2
    u = varargin{1};
    n = varargin{2};
    if round(n) ~= n || (n-1)/2 ~= floor(n/2)
        error('n must be an odd integer.')
    end
else
    error('Wrong number of inputs.')
end


% Displacement gradient
d = 1/12*n^2*(n^2-1);
nhalf = (n-1)/2;

[M, N] = size(u);
[dx, dy] = meshgrid(-(-nhalf:nhalf)/d,...
    -(-nhalf:nhalf)/d); % The extra - sign is due to a - sign appearing in the convolution equation
dudx = ones(M, N)*nan;
dudy = dudx;

dudx(1+nhalf:M-nhalf, 1+nhalf:N-nhalf) = convn(u, dx,'valid');
dudy(1+nhalf:M-nhalf, 1+nhalf:N-nhalf) = convn(u, dy,'valid');

