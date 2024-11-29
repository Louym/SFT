function [g,G]=genwindow(n,w)
i=1:w;
w1=gausswin(w,10)';
w2=2*sinc((i-w/2)/2048);
g=w1.*w2/352.512;
G=fft(g);
g1=[zeros(1,(n-w)/2) g zeros(1,(n-w)/2)];
G1=fft(g1);
figure
plot(w1)
title('时域窗函数');
xlabel('n');
ylabel('w[n]');
figure
plot(1./n.*(-n/2:n/2-1),abs(fftshift(G1)),"LineWidth",0.5)
title('频域窗函数');
xlabel('归一化频率');
ylabel('幅频响应');
end