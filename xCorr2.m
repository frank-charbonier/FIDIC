function C = xCorr2(A,B)
%xCorr2 FFT-based cross correlation
% 
% C = xCorr2(A,B)
% A and B are intput arrays; C is cross correlation.

A = fft2(A);
B = fft2(B);
B = conj(B);
C = A.*B;
C = ifft2(C);
C = real(C);
C = fftshift(C);
