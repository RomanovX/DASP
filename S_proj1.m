clc;
clear;
[s,fs]=wavread('F:\TU Delft courses\Speech and Audio Processing\clean.wav');
n1=wavread('F:\TU Delft courses\Speech and Audio Processing\noise1.wav');
spf=512;
y1=s+n1;
% Appending the speech signal with zeros
m=length(y1);
m1=m;
while(rem(m1,spf)~=0)
    m1=m1+1;
end
a=zeros(m1,1);
a(1:m,1)=y1;
% Segmentation and computing the 512 point DFT
f=m1/spf;
y=zeros(spf,f);
for i=1:f
    y(1:spf,i)=a(spf*(i-1)+1:spf*i,1);
end
Y=fft(y);
% Computing the periodogram of the individual segments
Py=(abs(Y).*abs(Y))/(spf);
% Computing the Bartlett estimate of the signal
PyB=sum(Py,2);
PyB=PyB/f;
plot(PyB)


        


        
