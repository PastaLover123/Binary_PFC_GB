
function [strain] = strainFilter(k,k2,G,sigma_0,sigma_1,gamma_1,n)
% strainFilter calculates a representation of the local hydrostatic strain
% for a given reference wave-vector given an input density field. This is
% done by identifying local shifts away from the reference wave-vectors by
% multiplying the Fourier transform of the density field with a Gaussian
% in Fourier space with a with sigma_1, centered about the reference
% wave-vector. Positive and negative deviations are then given a sign shift
% by multiplying the output signal by a linear function centered about the
% reference wave-vector. This signal is then convolved with the density
% again to identify regions where it is prominent and smooth with a
% low-pass filter with width sigma_0. The final signal amplitude is then
% modified by gamma_1.

strain = gamma_1*ifft2(k.*exp(-k2/sigma_0).*fft2(n.*(ifft2((k-G).*exp(-(k-G).^2/sigma_1).*fft2(n) )) ) );
end
