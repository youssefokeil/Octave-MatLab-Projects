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
fc=200; #10KHz

%%%%%% Time and frequency definition %%%%%%
# time step
ts= 0.1/fc;
# T is the simulation time
T=100;
# Number of time samples
N=ceil(T/ts);
# time vector
t= 0:ts:((N-1)*ts);

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

%%%%%% Baseband signal definition in Time %%%%%%
# m(t) definition
m= besselj(0,t).*cos(2*pi*10*t);
m(t>10)=0;
# plotting m(t) with t
figure(1);
plot(t,m);
xlabel("Time (sec)");
ylabel("m(t)");
box off;

%%%%%% frequency components of baseband sigmal %%%%%%
# fast fourier transform of baseband signal, scaled by 1/N for periodic signals, *ts for non-periodic signals
M=fftshift(fft(m)*ts);
# plotting frequency component of signal
figure(2);
plot(f, abs(M));
xlabel("Frequency (Hz)");
ylabel("M(f)");
box off;

%%%%%%%%% Calculating Total Energy %%%%%%%%
Energy_from_freq=sum(abs(M).^2)*df;


%%%%%% Finding Band Width %%%%%%%%%%
# we first have to find f=0
index_f0=find(f==0);
# Initializing accumulate energy
Energy_acc=0;
for index_f= index_f0:length(f)
  # we multiply by 2 to use the property of even function
  Energy_acc+=2*(abs(M(index_f)).^2)*df;
  # comparing it to 99% of total energy
  if(Energy_acc>=0.99*Energy_from_freq)
    BW=f(index_f);
    break
  end
end

%%%%%% DSB-SC Modulated signal definition in time %%%%%%
# s(t)
c=cos(2*pi*fc*t);
s= m.*c;
# plotting modulated signal with respect to time
figure(3);
plot(t,s);
xlabel("Time (sec)");
ylabel("s(t)")
box off;

%%%%%% DSB-SC Modulated signal definition in frequency %%%%%%
# S(f) using fft
S=fftshift(fft(s)*ts);
# plotting modulated signal with respect to frequency
figure(4);
plot(f,abs(S));
xlabel("Frequency (Hz)");
ylabel("S(f)")
box off;

%%%%%%%%% Band Pass Filter %%%%%%%%%%
H= zeros(size(f));
H(f>fc & f<fc+BW)=1;
H(f<-fc & f>-fc-BW)=1;
# plotting filter with respect to frequency
figure(5);
plot(f,abs(H));
xlabel("Frequency (Hz)");
ylabel("|H(f)|")
box off;

%%%%%% SSB (USB is transmitted) Modulated signal definition in frequency %%%%%%
# SSB frequency domain
SSB=S.*H;
# plotting modulated signal with respect to frequency
figure(6);
plot(f,abs(SSB));
xlabel("Frequency (Hz)");
ylabel("SSB |S(f)|")
box off;

%%%%%%%%%%% SSB signal time domain %%%%%%%%%%%
s=real(ifft(ifftshift(S))/ts);

%%%%%% Coherent detector (receiver) time response %%%%%%
# r(t)
g_received=s.*c;

%%%%%%%%%%%% Freq response %%%%%%%%%%%%%%%
G_received=fftshift(fft(g_received)*ts);

# plotting receiver signal with respect to time
figure(7);
plot(f,abs(G_received));
xlabel("Frequency (Hz)");
ylabel("G(f)");
box off;

%%%%%% Filter (ideal LPF) frequency response %%%%%%
# H(f)
H=abs(f)<BW;

%%%%%%%%%% Message after filter %%%%%%%%%%%%%
G_final=H.*G_received;
g_final=real(ifft(ifftshift(G_final)/ts));
# normalizing the curve
g_norm=g_final/max(g_final);

%%%%%%%% plotted signal after LPF %%%%%%%%%%
figure(8);
plot(t,g_norm,"--b");
xlabel("Time (sec)");
ylabel("g(t)")
box off;

%%%%%%% Plot Original Signal %%%%%%%%%%%%%%%
hold on;
plot(t,m/max(m),"-r");
legend("Signal after LPF","Original Signal");

