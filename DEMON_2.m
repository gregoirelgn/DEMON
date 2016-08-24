clear all;
close all;
clc;

%% loading data

[y,Fs]=audioread('2011_04_26_15_34_22_######_Ch0_00002.wav');

%% square method for extracting the envellop

%  bandpass filter
wc=5050; % wc is the center of the passing band (empirically chosen)

l=10e3; % l is the width of the bandpass filter
% Put y in a line vector
y=transpose(y);
% cutoff frequencies of the filter
Wc=[(wc-l/2) (wc+l/2)];
[b,a] = butter(1,Wc./Fs,'bandpass'); % Butterworth's filter, 1 order 
% filtering step
signal_mod=filter(b,a,y);

% Squaring operation

s=signal_mod.*signal_mod;

% low pass filter
W=500; % cutoff frequency of the law pass filter (empirically chosen)
[b,a] = butter(4,2.*W./Fs,'low'); % Butterworth's filter, 4 order 
s_filter=filter(b,a,s); % filter the squared signal

% square root of the filtered signal
env=sqrt(s_filter);



%% Diminution of the sampling frequency
% p is the rate of decreasing the sampling frequency
p=20;
n=length(env);
signal=zeros(1,floor(n/p));

for k=[1:floor(n/p)] 
    signal(k)=env(p*k);
end

fs=Fs/p; % fs is the new samplimg frequency

%% Ploting part

% plot the DEMON spectrum
subplot(211);
X=fftshift(fft(signal));
freq= linspace(0,Fs/2,length(signal)/2);
plot(freq, abs(X(length(signal)/2+1:length(signal))));
title('DEMON spectrum')
xlim([1 250]);
ylabel('|X(f)|');
xlabel('frequency');

% ploting the DEMON spectrogramm
subplot(212);
nfft=8192;
[S,F,T]=spectrogram(signal,hamming(nfft),nfft/2,nfft,fs);
pcolor(F,T,20*log10(abs(transpose(S)+3))); % waterfall plot % we add 3 to the spectrogram in order to improve the contrast  
title(['nfft = ',num2str(nfft)]);
xlim([0 250]); % we are only intersested in low frequencies
ylabel('time');
xlabel('frequency (Hz)');
shading interp;
colorbar;
