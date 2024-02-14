
clear
format long
N = 1024;

sobol = net(sobolset(256), N);
vd(:,1) = vdcorput(N-1,2);
vd(:,2) = vdcorput(N-1,4);
vd(:,3) = vdcorput(N-1,8);
vd(:,4) = vdcorput(N-1,16);
vd(:,5) = vdcorput(N-1,32);
vd(:,6) = vdcorput(N-1,64);
vd(:,7) = vdcorput(N-1,128);
vd(:,8) = vdcorput(N-1,256);
vd(:,9) = vdcorput(N-1,512);
vd(:,10) = vdcorput(N-1,1024);

%z1 = sqrt(-2.*(log(vd(1:end,1)))).*sinc(2*pi.*vd(1:end,8));
% z1 = sqrt(-2.*(log(vd(1:end,1)))).*tanh(vd(1:end,1));
% z1 = z1 - floor(z1);

%[~,lfval] = LFSR3([true false true true true false false false false false],N/2,N); %N=1024
%[~,lfval] = LFSR3([true false true false true false false false false],N/2,N); %N=512
%LF3(LF3 == -1) = 0;
%lfval = lfval/N;
lfval = rand(1,N);

sinc_sob = zeros(1,N);
abs_sinc_sob = zeros(1,N);
sinc_vdc = zeros(1,N);
abs_sinc_vd = zeros(1,N);
sinc_lfsr = zeros(1,N);
abs_sinc_lf = zeros(1,N);

EE = zeros(1,N);

X2_stream_sobol = zeros(N, N);
X2_stream_vdc = zeros(N, N);
X3_stream_vdc = zeros(N, N);
X4_stream_vdc = zeros(N, N);
X5_stream_vdc = zeros(N, N);
X2_stream_lfsr = zeros(N, N);

for i = 1:N
    for k = 1:N
        %if i/N > sobol(k,4)
        if i/N > vd(k,3)%3
            X2_stream_sobol(i,k) = 1;
        end
        if i/N > vd(k,8)%8
        %if i/N > z1(k)
            X2_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,5)%5
        %if i/N > z1(k)
            X3_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,10)%10
        %if i/N > z1(k)
            X4_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,9)
        %if i/N > z1(k)
            X5_stream_vdc(i,k) = 1;
        end
        if i/N > lfval(k)
            X2_stream_lfsr(i,k) = 1;
        end
    end
end

%N=1024, 1/42->24 1/20->51 1/6->170
%N=512, 1/42->12 1/20->26 1/6->85
%N=256, 1/42->6 1/20->13 1/6->42

% n1_s = and(X2_stream_sobol(102,:), circshift(X2_stream_sobol(102,:),3));
% n2_s = not(and(n1_s, X2_stream_sobol(12,:)));
% n3_s = not(and3(X2_stream_sobol(26,:), n2_s, circshift(n1_s,1)));
% n4_s = not(and3(X2_stream_sobol(85,:), n3_s, circshift(n1_s,2)));
% y_s = and(n4_s, circshift(X2_stream_sobol(102,:),6));
% sinc_sobol = sum(y_s)/N
% 
% n1_v = and(X2_stream_vdc(102,:), circshift(X2_stream_vdc(102,:),3));
% n2_v = not(and(n1_v, X2_stream_vdc(12,:)));
% n3_v = not(and3(X2_stream_vdc(26,:), n2_v, circshift(n1_v,1)));
% n4_v = not(and3(X2_stream_vdc(85,:), n3_v, circshift(n1_v,2)));
% y_v = and(n4_v, circshift(X2_stream_vdc(102,:),6));
% sinc_vdc = sum(y_v)/N
% 
% n1_lf = and(X2_stream_lfsr(102,:), circshift(X2_stream_lfsr(102,:),3));
% n2_lf = not(and(n1_lf, X2_stream_lfsr(12,:)));
% n3_lf = not(and3(X2_stream_lfsr(26,:), n2_lf, circshift(n1_lf,1)));
% n4_lf = not(and3(X2_stream_lfsr(85,:), n3_lf, circshift(n1_lf,2)));
% y_lf = and(n4_lf, circshift(X2_stream_lfsr(102,:),6));
% sinc_lfsr = sum(y_lf)/N


for i = 1:N
    %EE(i) = sin((i-1)/N);
    EE(i) = 1 - ((i/N)^2)/factorial(3) + ((i/N)^4)/factorial(5) - ((i/N)^6)/factorial(7);

    n1_s(i,:) = and(X2_stream_sobol(i,:), circshift(X2_stream_sobol(i,:),3));
    n2_s(i,:) = not(and(n1_s(i,:), X2_stream_sobol(24,:)));
    n3_s(i,:) = not(and3(X2_stream_sobol(51,:), n2_s(i,:), circshift(n1_s(i,:),1)));
    n4_s(i,:) = not(and3(X2_stream_sobol(170,:), n3_s(i,:), circshift(n1_s(i,:),2)));
    y_s(i,:) = n4_s(i,:);
    sinc_sobol(i) = sum(y_s(i,:))/N;
    abs_sinc_sob(i) = abs(sinc_sobol(i) - EE(i));

    %n1_v(i,:) = and(X2_stream_vdc(i,:), circshift(X2_stream_vdc(i,:),3));
    %n1_v(i,:) = and(X2_stream_lfsr(i,:), circshift(X2_stream_lfsr(i,:),2));
    n1_v(i,:) = and(X2_stream_sobol(i,:), circshift(X2_stream_sobol(i,:),2));%best with 15 Delays
    %n1_v(i,:) = and(X4_stream_vdc(i,:), circshift(X4_stream_vdc(i,:),3));
    n2_v(i,:) = not(and(n1_v(i,:), X2_stream_vdc(24,:)));
    n3_v(i,:) = not(and3(X3_stream_vdc(51,:), n2_v(i,:), circshift(n1_v(i,:),0)));
    n4_v(i,:) = not(and3(X4_stream_vdc(170,:), n3_v(i,:), circshift(n1_v(i,:),0)));
    y_v(i,:) = n4_v(i,:);
    sinc_vdc(i) = sum(y_v(i,:))/N;
    abs_sinc_vd(i) = abs(sinc_vdc(i) - EE(i));

    n1_lf(i,:) = and(X2_stream_lfsr(i,:), circshift(X2_stream_lfsr(i,:),3));
    n2_lf(i,:) = not(and(n1_lf(i,:), circshift(X2_stream_lfsr(24,:),1)));
    n3_lf(i,:) = not(and3(circshift(X2_stream_lfsr(51,:),1), n2_lf(i,:), circshift(n1_lf(i,:),1)));
    n4_lf(i,:) = not(and3(circshift(X2_stream_lfsr(170,:),1), n3_lf(i,:), circshift(n1_lf(i,:),2)));
    y_lf(i,:) = n4_lf(i,:);
    sinc_lfsr(i) = sum(y_lf(i,:))/N;
    abs_sinc_lf(i) = abs(sinc_lfsr(i) - EE(i));
end
MAE_sobol = mean(abs_sinc_sob)
MAE_vdc = mean(abs_sinc_vd)
MAE_lfsr = mean(abs_sinc_lf)