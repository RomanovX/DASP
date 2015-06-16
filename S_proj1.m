clc;
clear;
[s,fs]=wavread('F:\TU Delft courses\Speech and Audio Processing\clean.wav');
n1=wavread('F:\TU Delft courses\Speech and Audio Processing\noise1.wav');
spf=400;
y1=s+n1;
[m,n]=size(y1);
f=round(m/spf);
y=zeros(512,f);
for i=1:f
    y(1:spf,i)=y1(spf*(i-1)+1:spf*i,1);
end
y(spf+1:spf+m-(spf*f),f)=y1((f*spf)+1:m,1);
