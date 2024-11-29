function [I1, signal_evalfft]=innerloop(signal,n,w,B,k,d,g,G)
% 内循环，信号、信号长度、w、B、稀疏度

% (1) Choose a random σ, τ ∈ [n] with σ odd.
tao=randi([0 (n-1)]);
sigma=2*randi([1 ceil(n/2-1/2)])-1;
% (2.1) 时域重排
signal_permutation=zeros(1,w);
for i=((n-w)/2):((n+w)/2-1)
    signal_permutation(i-(n-w)/2+1)=signal(mod(sigma*i+tao,n)+1);
end
% figure
% subplot(3,1,1)
% plot(abs(fft(signal_permutation)));
% (2.2) 加窗
signal_window=signal_permutation.*g;
% (3)  Compute ˆzi = ˆyi(n/B) for i ∈ [B] 频域降采样
z=zeros(1,B);
for i=0:(B-1)
    for j=0:(ceil(w/B)-1)
        z(i+1)=z(i+1)+signal_window(mod(i+j*B,n)+1);
    end
end
z_dft=fft(z);
% subplot(3,1,2)
% plot(abs(fft(signal_window)));
% subplot(3,1,3)
% plot(abs(z_dft));
% (4) Hash function and offset
hsigma=mod(round(sigma*(0:(n-1))*B/n),B);
% (5) Location loops
[~,J] = maxk(abs(z_dft),d*k);
J=J-1;
I1=[];

for i=1:k*d
    ff=find(hsigma==J(i))-1;
    if ~isempty(ff)
        I1=[I1 ff];
    end
end
%(6) Estimation loop
if w==n
    offset_dc=1200;
else
    offset_dc=1.2*10^(-13)*sqrt(n/2^10);
end
signal_evalfft=z_dft(hsigma(I1+1)+1).*exp(I1.*1i*tao*2*pi/n)./(G(mod((sigma.*I1-hsigma(I1+1).*n/B),w)+1))*offset_dc;
end