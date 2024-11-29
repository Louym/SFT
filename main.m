clc
clear
close all

ctime = datestr(now,30);
tseed = str2double(ctime((end-5):end))/10+randi(2000);
rng(tseed)

n=2^22;
k=100;

%产生信号
frequency=rand(1,k/2)/2;
t=1:n;
signal=zeros(1,n);
for i=1:k/2
    signal=signal+(randi(10)+5)*cos(frequency(i)*2*pi.*(t-1))+(randi(10)+5)*sin(frequency(i)*2*pi.*(t-1));
end
figure
plot(0:(n-1),signal,"LineWidth",0.5)
xlim([0,n-1]);
xlabel("n")
ylabel("x[n]")
title('原信号')



disp(['FFT实现结果']);
tic
signal_fft=fft(signal);
toc
%绘图
figure
subplot(2,1,1)
plot(1./n.*(-n/2:n/2-1),abs(fftshift(signal_fft)),"LineWidth",0.5)
xlabel("归一化频率")
ylabel("|fft(X)|")
title('FFT计算的幅度谱')
subplot(2,1,2)
plot(1./n.*(-n/2:n/2-1),angle(fftshift(signal_fft)),"LineWidth",0.5)
xlabel("归一化频率")
ylabel("arg(fft(X))")
title('FFT计算的相位谱')

disp(['sFFT实现结果']);
%设置参数
% epsilon=0.265;%归一化截止频率
% delta=0.01;%峰值波纹
L=log2(n);
d=2;
B=round(sqrt(n*k));%/epsilon/log2(n/delta)
B=2^round(log2(B));
w=B*log2(n);
if w>n
    w=n;
end
w=n;
%准备窗函数
[g,G]=genwindow(n,w);

tic
%开始计算
signal_sfft=outerloop(signal,n,L,B,k/2,d,g,G,w);
toc
%绘图
figure
subplot(2,1,1)
plot(1./n.*(-n/2:n/2-1),abs(fftshift(signal_sfft)),"LineWidth",0.5)
xlabel("归一化频率")
ylabel("|sft(X)|")
title('SFT计算的幅度谱')
subplot(2,1,2)
plot(1./n.*(-n/2:n/2-1),angle(fftshift(signal_sfft)),"LineWidth",0.5)
xlabel("归一化频率")
ylabel("arg(sft(X))")
title('SFT计算的相位谱')