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
fc=1e4; #10KHz

%%%%%% Time and frequency definition %%%%%%
# time step
ts= 0.001;
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
m= besselj(0,t);
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

%%%%% Parseval Theorem Check %%%%%%%
# checking if energy of signal calculated from frequency domain equals energy of signal caalculated from time domain
Energy_from_time=sum(abs(m).^2)*ts; #absolute to get rid of imaginary numbers
Energy_from_freq=sum(abs(M).^2)*df;
# by looking at both variables values in workspace we can verify parseval's theorem.

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

%%%%%% Finding Amplitude modulation %%%%%%%%%%
# we first have to find fc, we'll assume fc=10*BW
fc_am=10*BW;
# 1+ka*m(t)>0 for all times
Energy_acc=0;
# nottice that this is the maximum value of ka, it will let (1+ka*m)=0
ka_max=-1/min(m);
ka=ka_max;
%% we'd want to graph envelope to check (1+ka*m)>0
envelope=(1+ka*m);
figure(3);
plot(t,envelope);
xlabel("Time (sec)");
ylabel("Envelope Amplitude");
box off;

%% Seems it's working fine, let's make our AM modulation
am_carrier=cos(2*pi*fc_am*t);
am_modulated=envelope.*am_carrier;
envelope=(1+ka*m);
figure(4);
plot(t,am_modulated);
xlabel("Time (sec)");
ylabel("AM modulated signal s(t)");
box off;

%%%%%% frequency components of AM modulated sigmal %%%%%%
# fast fourier transform of AM modulated signal, scaled by 1/N for periodic signals, *ts for non-periodic signals
S_am_modulated=fftshift(fft(am_modulated)*ts);
# plotting frequency component of signal
figure(5);
plot(f, abs(S_am_modulated));
xlabel("Frequency (Hz)");
ylabel("S(f)");
box off;

%%%%%% demodulated signal after ideal envelope detector %%%%%%
# am_rec(t), using hilbert function (envelope detector)
am_rec=abs(hilbert(am_modulated));
# plotting received signal after filter with respect to time
figure(6);
plot(t,am_rec,"--r");
# adding original song compared to received one
hold on;
plot(t,m,"--b");
# adding envelope signal after ka
hold on;
plot(t, envelope, "-k");
legend("Received Signal","Original Signal","Envelope Signal (after Ka)");
xlabel("Time (sec)");
ylabel("m(t)")
box off;
%%%%%%%%%% By noticimng the plot, the hilbert envelope function has traced the transmitted message after ka %%%%%%%%%%%%%%%%%%
