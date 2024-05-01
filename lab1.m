#mathematical operations
###using casting
##summation=a+b
##summ2=a+c
##aa=a
##mult_elem=a.*aa
##div_elem=a.\aa;
##multiply=a*b;
##a2=a*2;

% signal
Ts=0.01;
N=1000;
t= 0:Ts:(N-1)*Ts;
#t=linspace(0,10,11)
y1=cos(2*pi*t);

xlabel("time")
ylabel("signal")
title("cosine waveform")
y2=sin(2*pi*t)

##figure(1)
##plot(t,x, "--b");
##
##figure(2);
##plot( t, y2, "-r")

% freq domain
figure(3)
fs=1/Ts;
df=fs/N;
freq=-fs/2:df:fs/2-df; % N is even
##freq=-fs/2+0.5*df:df:fs/2-0.5*df; % N is odd
y1_f=fft(y1);
y2_f=fftshift(fft(y2)/N);
plot(freq,abs(y2_f));

