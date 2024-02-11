
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

%z1 = sqrt(-2.*(log(vd(1:end,1)))).*sin(2*pi.*vd(1:end,8));
% z1 = sqrt(-2.*(log(vd(1:end,1)))).*tanh(vd(1:end,1));
% z1 = z1 - floor(z1);

%[~,lfval] = LFSR3([true false true true true false false false false false],N/2,N); %N=1024
%[~,lfval] = LFSR3([true false true false true false false false false],N/2,N); %N=512
%LF3(LF3 == -1) = 0;
%lfval = lfval/N;
lfval = rand(1,N);

tanh_sobol = zeros(1,N);
abs_tanh_sob = zeros(1,N);
tanh_vdc = zeros(1,N);
abs_tanh_vd = zeros(1,N);
tanh_lfsr = zeros(1,N);
abs_tanh_lf = zeros(1,N);

EE = zeros(1,N);

X2_stream_sobol = zeros(N, N);
X3_stream_sobol = zeros(N, N);
X2_stream_vdc = zeros(N, N);
X3_stream_vdc = zeros(N, N);
X4_stream_vdc = zeros(N, N);
X5_stream_vdc = zeros(N, N);
X2_stream_lfsr = zeros(N, N);

for i = 1:N
    for k = 1:N
        %if i/N > sobol(k,4)
        if i/N > vd(k,3)
            X2_stream_sobol(i,k) = 1;
        end
        if i/N > sobol(k,220)
        %if i/N > vd(k,10)
            X3_stream_sobol(i,k) = 1;
        end
        if i/N > vd(k,1)
        %if i/N > z1(k)
            X2_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,10)
        %if i/N > z1(k)
            X3_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,8)
        %if i/N > z1(k)
            X4_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,4)
        %if i/N > z1(k)
            X5_stream_vdc(i,k) = 1;
        end
        if i/N > lfval(k)
            X2_stream_lfsr(i,k) = 1;
        end
    end
end

%N=1024, 62/153->415 17/42->414 2/5->409 1/3->341
%N=512, 1/42->12 1/20->26 1/6->85
%N=256, 1/42->6 1/20->13 1/6->42

% n1_s = and(X2_stream_sobol(102,:), circshift(X2_stream_sobol(102,:),3));
% n2_s = not(and(n1_s, X2_stream_sobol(12,:)));
% n3_s = not(and3(X2_stream_sobol(26,:), n2_s, circshift(n1_s,1)));
% n4_s = not(and3(X2_stream_sobol(85,:), n3_s, circshift(n1_s,2)));
% y_s = and(n4_s, circshift(X2_stream_sobol(102,:),6));
% tanh_sobol = sum(y_s)/N
% 
% n1_v = and(X2_stream_vdc(102,:), circshift(X2_stream_vdc(102,:),3));
% n2_v = not(and(n1_v, X2_stream_vdc(12,:)));
% n3_v = not(and3(X2_stream_vdc(26,:), n2_v, circshift(n1_v,1)));
% n4_v = not(and3(X2_stream_vdc(85,:), n3_v, circshift(n1_v,2)));
% y_v = and(n4_v, circshift(X2_stream_vdc(102,:),6));
% tanh_vdc = sum(y_v)/N
% 
% n1_lf = and(X2_stream_lfsr(102,:), circshift(X2_stream_lfsr(102,:),3));
% n2_lf = not(and(n1_lf, X2_stream_lfsr(12,:)));
% n3_lf = not(and3(X2_stream_lfsr(26,:), n2_lf, circshift(n1_lf,1)));
% n4_lf = not(and3(X2_stream_lfsr(85,:), n3_lf, circshift(n1_lf,2)));
% y_lf = and(n4_lf, circshift(X2_stream_lfsr(102,:),6));
% tanh_lfsr = sum(y_lf)/N


for i = 1:N
    %EE(i) = sin((i-1)/N);
    EE(i) = (i/N) - ((i/N)^3)/3 + 2*((i/N)^5)/15 - 17*((i/N)^7)/315 + 62*((i/N)^9)/2835;

    n1_s(i,:) = and(X2_stream_sobol(i,:), circshift(X2_stream_sobol(i,:),4));
    n2_s(i,:) = not(and(n1_s(i,:), circshift(X2_stream_sobol(415,:),0)));
    n3_s(i,:) = not(and3(circshift(X3_stream_sobol(414,:),0), n2_s(i,:), circshift(n1_s(i,:),1)));
    n4_s(i,:) = not(and3(circshift(X2_stream_sobol(409,:),0), n3_s(i,:), circshift(n1_s(i,:),2)));
    n5_s(i,:) = not(and3(circshift(X3_stream_sobol(341,:),0), n4_s(i,:), circshift(n1_s(i,:),0)));
    y_s(i,:) = and(n5_s(i,:), circshift(X2_stream_sobol(i,:),0));
    tanh_sobol(i) = sum(y_s(i,:))/N;
    abs_tanh_sob(i) = abs(tanh_sobol(i) - EE(i));

    input = X2_stream_sobol(i,:);
    %input = X2_stream_vdc(i,:);
    %input = X2_stream_lfsr(i,:);
    n1_v(i,:) = and(input, circshift(input,3));
    n2_v(i,:) = not(and(n1_v(i,:), X2_stream_vdc(415,:)));
    n3_v(i,:) = not(and3(X3_stream_vdc(414,:), n2_v(i,:), circshift(n1_v(i,:),0)));
    n4_v(i,:) = not(and3(X4_stream_vdc(409,:), n3_v(i,:), circshift(n1_v(i,:),0)));
    n5_v(i,:) = not(and3(X5_stream_vdc(341,:), n4_v(i,:), circshift(n1_v(i,:),0)));
    y_v(i,:) = and(n5_v(i,:), circshift(input,0));
    tanh_vdc(i) = sum(y_v(i,:))/N;
    abs_tanh_vd(i) = abs(tanh_vdc(i) - EE(i));

    n1_lf(i,:) = and(X2_stream_lfsr(i,:), circshift(X2_stream_lfsr(i,:),4));
    n2_lf(i,:) = not(and(n1_lf(i,:), circshift(X2_stream_lfsr(415,:),1)));
    n3_lf(i,:) = not(and3(circshift(X2_stream_lfsr(414,:),1), n2_lf(i,:), circshift(n1_lf(i,:),1)));
    n4_lf(i,:) = not(and3(circshift(X2_stream_lfsr(409,:),1), n3_lf(i,:), circshift(n1_lf(i,:),2)));
    n5_lf(i,:) = not(and3(circshift(X2_stream_vdc(341,:),1), n4_lf(i,:), circshift(n1_lf(i,:),3)));
    y_lf(i,:) = and(n5_lf(i,:), circshift(X2_stream_lfsr(i,:),8));
    tanh_lfsr(i) = sum(y_lf(i,:))/N;
    abs_tanh_lf(i) = abs(tanh_lfsr(i) - EE(i));
end
MAE_sobol = mean(abs_tanh_sob)
MAE_vdc = mean(abs_tanh_vd)
MAE_lfsr = mean(abs_tanh_lf)