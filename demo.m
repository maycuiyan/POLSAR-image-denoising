%% Simulate a polarimetric SAR (PolSAR) image from an optical image

A = double(imread('./data/lena.tif')); % amplitude image
L = 8;                          % number of looks
mrows = size(A,1);              % number of rows
ncols = size(A,2);              % number of columns

% PolSAR is a multi-channel (6 in total) complex image
T11 = zeros(mrows,ncols);
T12 = zeros(mrows,ncols);
T13 = zeros(mrows,ncols);
T22 = zeros(mrows,ncols);
T23 = zeros(mrows,ncols);
T33 = zeros(mrows,ncols);
for i = 1:mrows
    for j = 1:ncols
        Sigma = diag([A(i,j,3),A(i,j,1),A(i,j,2)]);
        X = Sigma*(randn(3,L) + 1i * randn(3,L))/sqrt(2);
        T = X*X'/L;
        T11(i,j) = T(1,1);
        T12(i,j) = T(1,2);
        T13(i,j) = T(1,3);
        T22(i,j) = T(2,2);
        T23(i,j) = T(2,3);
        T33(i,j) = T(3,3);
    end
end

%% Visualize the noisy image

I_noisy(:,:,1)=sqrt(T22);
I_noisy(:,:,2)=sqrt(T33);
I_noisy(:,:,3)=sqrt(T11);
subplot(1,2,1), imshow(I_noisy/mean(I_noisy(:))/2);

%% Filter the image

KAPPA = 50;     % smoothing factor
[Tf11,Tf12,Tf13,Tf22,Tf23,Tf33] = polsar_filter(T11,T12,T13,T22,T23,T33, KAPPA);

%% Visualize the filtered image

I_filtered(:,:,1)=sqrt(Tf22);
I_filtered(:,:,2)=sqrt(Tf33);
I_filtered(:,:,3)=sqrt(Tf11);
subplot(1,2,2), imshow(I_filtered/mean(I_filtered(:))/2);
