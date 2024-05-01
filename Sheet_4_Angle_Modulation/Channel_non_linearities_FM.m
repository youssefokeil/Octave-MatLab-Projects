
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

%%%%%%%% FM Modulated Signal %%%%%%%%%%%%%
# time domain, cumtrapz is integration function
s= Ac*cos(2*pi*fc*t+2*pi*kf*cumtrapz(t,m));
# frequency domain
S=fftshift(fft(s)/N);
%%%%%%%%%%% plotting FM modulated signal %%%%%%%%
# plotting modulated signal with respect to time
figure(1);
plot(t,s);
xlabel("Time (sec)");
ylabel("s(t)")
# plotting modulated signal with respect to frequency
figure(2);
# using stem instead of plot
stem(f,abs(S));
xlabel("Frequency (Hz)");
ylabel("|S(f)|")
box off;


%%%%%%%%%% Channel Non-linearities %%%%%%%%%%
s_out= s + 0.5 * s.^2 + s.^3;
S=fftshift(fft(s)/N);
S_out=fftshift(fft(s_out)/N);

# plotting modulated signal with respect to frequency
figure(3);
# using stem instead of plot
stem(f,abs(S_out));
xlabel("Frequency (Hz)");
ylabel("|S(f)| after Non-Linear Channel")
box off;

%%%%%% using carson's law to get Band Width %%%%%%%%%
df=kf*max(m);
BW=2*fm+2*df;


%%%%%%%% FM receiver %%%%%%%%%%%%%%%%%%%
%%% Band Pass Filter %%%%%%%%%
H=zeros(size(f));
H(f>(fc-0.5*BW)&f<(fc+0.5*BW))=1;
H(f<-(fc-0.5*BW)&f>-(fc+0.5*BW))=1;

%%%%% Retrieve FM signal %%%%%%%%
S_out=H.*S_out;
# plotting modulated signal with respect to frequency
figure(4);
# using stem instead of plot
stem(f,abs(S));
xlabel("Frequency (Hz)");
ylabel("|S(f)| after Band Pass Filter")
box off;
% plotting signal time domain after retrieval %%
s_out=real(ifft(ifftshift(S_out)*ts));

%%%% Plotting original signal %%%%%%%
%%%%%% calling figure 1 again which had the modulated signal with respect tto time
figure(1)
box off;
hold on;
%%%%% Plotting retrieved signal %%%%%%%%%%
plot(t,s_out/max(s_out),"--r");
legend("Original Signal", "Retrieved Signal");


