# clear command window
clc;
# remove all variables from workspace
clear
# close all figures
close all

%%%%%% Frequencies initialization %%%%%%
# carrier frequency initialization
fc=1e4; #10KHz
# message frequency initialization
fm=100; #100Hz
Am=10;   %try Am=10
Ac=1;
kf=100;  %try kf=100
beta=kf*Am/fm;

%%%%%% Time and frequency definition %%%%%%
# Sampling time is identified using fc so that we can sample carrier signal tooo
ts= 0.1/fc;
# T is the last value of t on graph, we're drawing 2 periods of the message signal
T=2/fm;
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
  f= -(0.5*fs): df : (0.5*fs-df)
else #odd
  f= -(0.5*fs-0.5*df): df : (0.5*fs-0.5*df)
end

%%%%%% Baseband signal definition in Time %%%%%%
# m(t) definition
m= Am*sin(2*pi*fm*t);
# plotting m(t) with t
figure(1);
plot(t,m);
xlabel("Time (sec)");
ylabel("m(t)");
box off;

%%%%%% frequency components of baseband sigmal %%%%%%
# fast fourier transform of baseband signal
M=fftshift(fft(m)/N);
# plotting frequency component of signal
figure(2);
plot(f, abs(M));
xlabel("Frequency (Hz)");
ylabel("|M(f)|");
box off;

%%%%%%%% Modulated Signal %%%%%%%%%%%%%
# time domain, cumtrapz is integration function
s= Ac*cos(2*pi*fc*t+2*pi*kf*cumtrapz(t,m));
# frequency domain
S=fftshift(fft(s)/N);

%%%%%%%%%%% plotting FM modulated signal %%%%%%%%
# plotting modulated signal with respect to time
figure(3);
plot(t,s);
xlabel("Time (sec)");
ylabel("s(t)")
# plotting modulated signal with respect to frequency
figure(4);
# using stem instead of plot
stem(f,abs(S));
xlabel("Frequency (Hz)");
ylabel("|S(f)|")
box off;

%%%%%%%%%%%% Spectrum using formula of WideBand FM %%%%%%%%%%%%%%
n=-100:100;
S_analytical = 0.5*besselj(n,beta);
f_analytical = fc + n*fm;
hold all;
stem(f_analytical,abs(S_analytical),"--r");
legend("Simulation","Analytical");

%%%%%%%%% Calculating BW %%%%%%%%%%%%%%%%
%%%%% Using Inuversal Law %%%%%%%%
% BW where no component outside has a value more than 1% of carrier amplitude 𝑨𝒄 %
% index of fc carrier
index_fc=find((abs(f-fc)) == min(abs(f-fc)));
for index_f = length(f):-1:index_fc;
  if(abs(S(index_f))>=0.5*0.01*Ac)
    BW=2*abs(f(index_f)-fc);
    break;
  end
end
%%%%% Using Carson Law %%%%%%%%%%
% Carson’s rule (Empirical relation) %
df=kf*Am;
BW_carson=2*fm+2*df;
