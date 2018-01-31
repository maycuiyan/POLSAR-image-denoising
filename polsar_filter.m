function [Tf11,Tf12,Tf13,Tf22,Tf23,Tf33] = polsar_filter(T11,T12,T13,T22,T23,T33, KAPPA, ENL, R)

% This program implements polarimetric SAR (POLSAR) data filtering according to the following paper:
% Y. Cui, Y. Yamaguchi, H. Kobayashi, and J. Yang, "FILTERING OF POLARIMETRIC SYNTHETIC APERTURA RADAR IMAGES: A SEQUENTIAL APPROACH"
% presented at IGARSS2012, Munich.

% Argument Description:
% T11,T12,T13,T22,T23,T33: elements of POLSAR coherency/covariance matrix
% KAPPA: user-defined smoothing constant
% ENL: equivalent number of looks of the POLSAR data
% R: half size of the window for local ENL estimation


if nargin<6
    error('please input all the channels of the POLSAR coherency/covariance matrix data');
elseif nargin==6
    error('please input the user-defined smoothing constant');
elseif nargin==7
    l_enl=fPOLSARENL(T11,T12,T13,T22,T23,T33,3);
    [n,xout] = hist(l_enl(:),100);
    [maxCount, ind] = max(n);
    ENL = xout(ind);
    R=3;
elseif nargin==8
    R=3;
end
[Tf11,Tf12,Tf13,Tf22,Tf23,Tf33] = fPOLSARSMAP(double(T11),double(T12),double(T13),double(T22),double(T23),double(T33),double(ENL),uint8(R),double(KAPPA));