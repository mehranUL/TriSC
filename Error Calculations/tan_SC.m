
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

tan_sobol = zeros(1,N);
abs_tan_sob = zeros(1,N);
tan_vdc = zeros(1,N);
cot_vdc = zeros(1,N);
abs_tan_vd = zeros(1,N);
tan_lfsr = zeros(1,N);
abs_tan_lf = zeros(1,N);

sin_sob = zeros(1,N);
abs_sin_sob = zeros(1,N);
sin_vdc = zeros(1,N);
abs_sin_vd = zeros(1,N);
sin_lfsr = zeros(1,N);
abs_sin_lf = zeros(1,N);

cos_sobol = zeros(1,N);
abs_cos_sob = zeros(1,N);
cos_vdc = zeros(1,N);
abs_cos_vd = zeros(1,N);
cos_lfsr = zeros(1,N);
abs_cos_lf = zeros(1,N);

EE = zeros(1,N);

X2_stream_sobol = zeros(N, N);
X3_stream_sobol = zeros(N, N);
X2_stream_vdc_cos = zeros(N, N);
X3_stream_vdc_cos = zeros(N, N);
X4_stream_vdc_cos = zeros(N, N);
X5_stream_vdc_cos = zeros(N, N);
X2_stream_vdc_sin = zeros(N, N);
X3_stream_vdc_sin = zeros(N, N);
X4_stream_vdc_sin = zeros(N, N);
X5_stream_vdc_sin = zeros(N, N);
X2_stream_lfsr = zeros(N, N);
X2_stream_vdc_sinc = zeros(N, N);
X3_stream_vdc_sinc = zeros(N, N);
X4_stream_vdc_sinc = zeros(N, N);

for i = 1:N
    for k = 1:N
        %if i/N > sobol(k,1)
        if i/N > vd(k,3)
            X2_stream_sobol(i,k) = 1;
        end
        if i/N > sobol(k,9)
        %if i/N > vd(k,10)
            X3_stream_sobol(i,k) = 1;
        end
        if i/N > vd(k,4) %3,6,7 %2,4
        %if i/N > z1(k)
            X2_stream_vdc_cos(i,k) = 1;
        end
        if i/N > vd(k,9)%3,9
        %if i/N > z1(k)
            X3_stream_vdc_cos(i,k) = 1;
        end
        if i/N > vd(k,4)%4
        %if i/N > z1(k)
            X4_stream_vdc_cos(i,k) = 1;
        end
        if i/N > vd(k,8)%8
        %if i/N > z1(k)
            X5_stream_vdc_cos(i,k) = 1;
        end
        if i/N > vd(k,6)%6
        %if i/N > z1(k)
            X2_stream_vdc_sin(i,k) = 1;
        end
        if i/N > lfval(k)
            X2_stream_lfsr(i,k) = 1;
        end
    end
end

%N=1024, 1/56->18 1/30->34 1/12->85 1/2->512
%N=512, 1/42->12 1/20->26 1/6->85
%N=256, 1/42->6 1/20->13 1/6->42
jj = 0;
ii = 1;

for i = 1:N
    %EE(i) = cos((i-1)/N);
    EE(i) = (i/N) + ((i/N)^3)/3 + (2/15)*((i/N)^5) + (17/315)*((i/N)^7);

%     n1_s(i,:) = and(X2_stream_sobol(i,:), circshift(X2_stream_sobol(i,:),4));
%     n2_s(i,:) = not(and(n1_s(i,:), circshift(X2_stream_sobol(18,:),0)));
%     n3_s(i,:) = not(and3(circshift(X3_stream_sobol(34,:),0), n2_s(i,:), circshift(n1_s(i,:),1)));
%     n4_s(i,:) = not(and3(circshift(X2_stream_sobol(85,:),0), n3_s(i,:), circshift(n1_s(i,:),2)));
%     y_s(i,:) = not(and3(circshift(X3_stream_sobol(512,:),0), n4_s(i,:), circshift(n1_s(i,:),3)));
%     tan_sobol(i) = sum(y_s(i,:))/N;
%     abs_tan_sob(i) = abs(tan_sobol(i) - EE(i));

    input = X2_stream_sobol(i,:);
    %input = X2_stream_vdc(i,:);
    %input = X2_stream_lfsr(i,:);
    n1_v_sin(i,:) = and(input, circshift(input,2));%best with 15 Delays
    %n1_v(i,:) = and(X4_stream_vdc(i,:), circshift(X4_stream_vdc(i,:),3));
    n2_v_sin(i,:) = not(and(n1_v_sin(i,:), X2_stream_vdc_sin(24,:)));
    n3_v_sin(i,:) = not(and3(X2_stream_vdc_sin(51,:), n2_v_sin(i,:), circshift(n1_v_sin(i,:),0)));
    n4_v_sin(i,:) = not(and3(X2_stream_vdc_sin(170,:), n3_v_sin(i,:), circshift(n1_v_sin(i,:),0)));
    y_v_sin(i,:) = and(n4_v_sin(i,:), circshift(input,0));
    sin_vdc(i) = sum(y_v_sin(i,:))/N;


    input = X2_stream_sobol(i,:);
    %input = X2_stream_vdc(i,:);
    %input = X2_stream_lfsr(i,:);
    n1_v_cos(i,:) = and(input, circshift(input,2));%best with 15 Delays
    n2_v_cos(i,:) = not(and(n1_v_cos(i,:), X2_stream_vdc_cos(18,:)));
    n3_v_cos(i,:) = not(and3(X3_stream_vdc_cos(34,:), n2_v_cos(i,:), circshift(n1_v_cos(i,:),0)));
    n4_v_cos(i,:) = not(and3(X4_stream_vdc_cos(85,:), n3_v_cos(i,:), circshift(n1_v_cos(i,:),0)));
    y_v_cos(i,:) = not(and3(X5_stream_vdc_cos(512,:), n4_v_cos(i,:), circshift(n1_v_cos(i,:),0)));
    cos_vdc(i) = sum(y_v_cos(i,:))/N;

    

    if sin_vdc(i) < cos_vdc(i)
        [~,y_v_tan(i,:)] = CORLD_DIV(N*sin_vdc(i), N*cos_vdc(i), N);
        tan_vdc(i) = sum(y_v_tan(i,:))/N;
        jj = jj+1;
        %abs_tan_vd(i) = abs(tan_vdc(i) - EE(i));
    else
        tmp = sin_vdc(i) - cos_vdc(i);
        while tmp > cos_vdc(i)
            tmp = sin_vdc(i) - cos_vdc(i);
            ii = ii+1;
        end
        jj = jj+1;
        [~,y_v_tan(jj,:)] = CORLD_DIV(N*tmp, N*cos_vdc(jj), N);
        tan_vdc(jj) = ii + sum(y_v_tan(i,:))/N;
        
    end
    abs_tan_vd(i) = abs(tan_vdc(i) - EE(i));
%     for xx = 1:jj
%         tan_vdc(jj) = 
%     end
   

%     n1_lf(i,:) = and(X2_stream_lfsr(i,:), circshift(X2_stream_lfsr(i,:),4));
%     n2_lf(i,:) = not(and(n1_lf(i,:), circshift(X2_stream_lfsr(18,:),1)));
%     n3_lf(i,:) = not(and3(circshift(X2_stream_lfsr(34,:),1), n2_lf(i,:), circshift(n1_lf(i,:),1)));
%     n4_lf(i,:) = not(and3(circshift(X2_stream_lfsr(85,:),1), n3_lf(i,:), circshift(n1_lf(i,:),2)));
%     y_lf(i,:) = not(and3(circshift(X2_stream_lfsr(512,:),1), n4_lf(i,:), circshift(n1_lf(i,:),3)));
%     tan_lfsr(i) = sum(y_lf(i,:))/N;
%     abs_tan_lf(i) = abs(tan_lfsr(i) - EE(i));
end
%MAE_sobol = mean(abs_tan_sob)
MAE_vdc = mean(abs_tan_vd)
%MAE_lfsr = mean(abs_tan_lf)