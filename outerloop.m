function x=outerloop(signal,n,L,B,k,d,g,G,w)
x=zeros(1,n);
y=zeros(1,n);
I=[];
signal_evalfft=[];
for i=1:L
    [I_temp,signal_evalfft_temp]=innerloop(signal,n,w,B,k,d,g,G);
    I=[I I_temp];
    signal_evalfft=[signal_evalfft,signal_evalfft_temp];
%     y(I_temp)=signal_evalfft;
%     figure
%     plot(abs(y));
end
table=tabulate(I);
[l,~]=size(table);
for i=1:l
    if table(i,2)>=L/2
        temp1=I==table(i,1);    
        temp2=signal_evalfft(temp1);
        x(table(i,1)+1)=median(real(temp2))+1i*median(imag(temp2));        
    end
end
end