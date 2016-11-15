function wx = gen_noise( SNR_dB, tx)

SNR = 10^(SNR_dB/10); % calculate linear SNR

Pa = sum(abs(tx).^2)/length(tx) % calculate signal power

L=length(tx);
wx = sqrt(0.5)*(randn(L, 1)+ 1i*randn(L, 1)); 
% computed noise, randn generated PDF with variance 1
% noise power is variance, so it is 0.5+0.5 now 
% (add noise power of real and complex part)

% SNR = Pa / Pw
% Desired noise power is Pa / SNR
% Noise has a quadratic relation with amplitude
% So scale with sqrt( desired_noise_power / noise power) 
% because urrent noise power is 1
wx = sqrt(Pa / SNR)*wx;

end