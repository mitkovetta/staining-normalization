function [Inorm, H, E] = normalizeStaining(I, Io, beta, alpha, HERef, maxCRef)
% Normalize the staining appearance of images originating from H&E stained
% sections.
%
% Example use:
% See test.m.
%
% Input:
% I         - RGB input image;
% Io        - (optional) transmitted light intensity (default: 240);
% beta      - (optional) OD threshold for transparent pixels (default: 0.15);
% alpha     - (optional) tolerance for the pseudo-min and pseudo-max (default: 1);
% HERef     - (optional) reference H&E OD matrix (default value is defined);
% maxCRef   - (optional) reference maximum stain concentrations for H&E (default value is defined);
%
% Output:
% Inorm     - normalized image;
% H         - (optional) hematoxylin image;
% E         - (optional)eosin image;
%
% References:
% A method for normalizing histology slides for quantitative analysis. M.
% Macenko et al., ISBI 2009
%
% Efficient nucleus detector in histopathology images. J.P. Vink et al., J
% Microscopy, 2013
%
% Copyright (c) 2013, Mitko Veta
% Image Sciences Institute
% University Medical Center
% Utrecht, The Netherlands
% 
% See the license.txt file for copying permission.
%

% transmitted light intensity
if ~exist('Io', 'var') || isempty(Io)
    Io = 240;
end

% OD threshold for transparent pixels
if ~exist('beta', 'var') || isempty(beta)
    beta = 0.15;
end

% tolerance for the pseudo-min and pseudo-max
if ~exist('alpha', 'var') || isempty(alpha)
    alpha = 1;
end

% reference H&E OD matrix
if ~exist('HERef', 'var') || isempty(HERef)
    HERef = [
        0.5626    0.2159
        0.7201    0.8012
        0.4062    0.5581
        ];
end

% reference maximum stain concentrations for H&E
if ~exist('maxCRef)', 'var') || isempty(maxCRef)
    maxCRef = [
        1.9705
        1.0308
        ];
end

h = size(I,1);
w = size(I,2);

I = double(I);

I = reshape(I, [], 3);

% calculate optical density
OD = -log((I+1)/Io);

% remove transparent pixels
ODhat = OD(~any(OD < beta, 2), :);

% calculate eigenvectors
[V, ~] = eig(cov(ODhat));

% project on the plane spanned by the eigenvectors corresponding to the two
% largest eigenvalues
That = ODhat*V(:,2:3);

% find the min and max vectors and project back to OD space
phi = atan2(That(:,2), That(:,1));

minPhi = prctile(phi, alpha);
maxPhi = prctile(phi, 100-alpha);

vMin = V(:,2:3)*[cos(minPhi); sin(minPhi)];
vMax = V(:,2:3)*[cos(maxPhi); sin(maxPhi)];

% a heuristic to make the vector corresponding to hematoxylin first and the
% one corresponding to eosin second
if vMin(1) > vMax(1)
    HE = [vMin vMax];
else
    HE = [vMax vMin];
end

% rows correspond to channels (RGB), columns to OD values
Y = reshape(OD, [], 3)';

% determine concentrations of the individual stains
C = HE \ Y;

% normalize stain concentrations
maxC = prctile(C, 99, 2);

C = bsxfun(@rdivide, C, maxC);
C = bsxfun(@times, C, maxCRef);

% recreate the image using reference mixing matrix
Inorm = Io*exp(-HERef * C);
Inorm = reshape(Inorm', h, w, 3);
Inorm = uint8(Inorm);

if nargout > 1
    H = Io*exp(-HERef(:,1) * C(1,:));
    H = reshape(H', h, w, 3);
    H = uint8(H);
end

if nargout > 2
    E = Io*exp(-HERef(:,2) * C(2,:));
    E = reshape(E', h, w, 3);
    E = uint8(E);
end

end
