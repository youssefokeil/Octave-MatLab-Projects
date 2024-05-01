
# clear command window
clc;
# remove all variables from workspace
clear
# close all figures
close all
# importing signal package
pkg load signal

%%%%%% Frequencies initialization %%%%%%
# carrier frequency initialization
fc=1e3; #1KHz
# message frequency initialization
fm=100; #100Hz
Am=1;
Ac=1;
kf=10;
ka=0.5;
beta=kf*Am/fm;

%%%%%% Time and frequency definition %%%%%%
# Sampling time is identified using fc so that we can sample carrier signal tooo
ts= 0.1/fc/10;
# T is the last value of t on graph, we're drawing 2 periods of the message signal
T=10/fm;
# Number of time samples
N=ceil(T/ts);
# time vector
t= 0:ts:((N-1)*ts)

%%%%%% Frequency domain %%%%%%
# df=fs/N
df=1/T;
# sampling frequency, inverse of sampling tinme
fs=1/ts;
# frequency vector
if(rem(N,2)==0) #even
  f= -(0.5*fs): df : (0.5*fs-df);
else #odd
  f= -(0.5*fs-0.5*df): df : (0.5*fs-0.5*df);
end



%%%%% Baseband Signal %%%%%%%
m= Am * cos(2*pi*fm*t);




%%%%%%%% AM Modulation %%%%%%%%%%
%%%%%% Finding Amplitude modulation %%%%%%%%%%

##### Am Modulation #######
s_am=Ac*(1+ka*m).*cos(2*pi*fc*t);
figure(5);
plot(t,s_am);
xlabel("Time (sec)");
ylabel("Amplitude Modulation s(t)");
box off;

%%%% AM frequency response %%%%%%%
figure(6);
S_am=fftshift(fft(s_am)/N);
stem(f,abs(S_am));
xlabel("Frequency (Hz)");
ylabel("Amplitude Modulation S(f)");
box off;

%%%%%%%%%% Channel Non-linearities %%%%%%%%%%
s_out= s_am + 0.5 * s_am.^2 + s_am.^3;
S_out=fftshift(fft(s_out)/N);
%%%% AM frequency response after channel non-linearities%%%%%%%
figure(7);
stem(f,abs(S_out));
xlabel("Frequency (Hz)");
ylabel("Amplitude Modulation S(f)");
box off;

%%%%%%% identifying BW of AM %%%%%%%%%%
BW=2*fm;
%%% Band Pass Filter %%%%%%%%%
H=zeros(size(f));
H(f>(fc-0.5*BW-df)&f<(fc+0.5*BW+df))=1;
H(f<-(fc-0.5*BW-df)&f>-(fc+0.5*BW+df))=1;
%%%%% Retrieval using band pass %%%%%%%%%%%
S_out=H.*S_out;

# plotting modulated signal with respect to frequency
figure(8);
# using stem instead of plot
stem(f,abs(S_out));
xlabel("Frequency (Hz)");
ylabel("|S(f)| after Band Pass Filter")
box off;
% plotting signal time domain after retrieval %%
s_out=real(ifft(ifftshift(S_out)*ts));
%%%%%% calling figure 5 again which had the modulated signal with respect tto time
figure(5)
hold on;
%%%%% Plotting retrieved signal %%%%%%%%%%
plot(t,s_out/max(s_out),"--r");
legend("Original Signal", "Retrieved Signal");
