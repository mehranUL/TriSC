%Calculating the sin(x) using SC calculations.
%The input range is x is in [0,1] interval (in Radian)

function [sin_vdc,sin_sobol,sin_lfsr] = Sine_vdc_chu(X, N)
format long
if X > 1 || X < 0
    fprintf("Error");
else
%clear

%N = 1024;
sobol = net(sobolset(256), N);

[~,lfval] = LFSR_TrigonoSC([true false true true true false false false false false],X,N); %N=1024
%[~,lfval] = LFSR_TrigonoSC([true false true false true false false false false],X,N); %N=512
lfval = lfval/N;

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

%sin_vdc = zeros(1,N);

X2_stream_vdc = zeros(1, N);
X2_stream_sobol = zeros(1, N);
X2_stream_lfsr = zeros(1, N);

X2_stream_vdc5040 = zeros(1, N);
X2_stream_vdc120 = zeros(1, N);
X2_stream_vdc6 = zeros(1, N);
%for i = 1:N
    for k = 1:N
        %if i/N > sobol(k,9)
        if X > vd(k,2)
            X2_stream_vdc(k) = 1;
        end
        if X > sobol(k,9)
            X2_stream_sobol(k) = 1;
        end
        if X > lfval(k)
            X2_stream_lfsr(k) = 1;
        end
        if 1/5040 > vd(k,7)
            X2_stream_vdc5040(k) = 1;
        end
        if 1/120 > vd(k,7)
            X2_stream_vdc120(k) = 1;
        end
        if 1/6 > vd(k,7)
            X2_stream_vdc6(k) = 1;
        end
    end
%end

    
    
    %input = X2_stream_sobol(i,:);
    input = X2_stream_vdc;
    n1_v = and(input, circshift(input,1));
    n2_v = not(and(n1_v, X2_stream_vdc5040));
    n3_v = and(n2_v, X2_stream_vdc120);
    n4_v = not(and(n3_v, circshift(n1_v, 0)));%1
    n5_v = and(n4_v, X2_stream_vdc6);
    n6_v = not(and(n5_v, circshift(n1_v, 0)));%2
    y_v = and(n6_v, circshift(input,0));%2
    sin_vdc = sum(y_v)/N;

    input = X2_stream_sobol;
    n1_s = and(input, circshift(input,1));
    n2_s = not(and(n1_s, X2_stream_vdc5040));
    n3_s = and(n2_s, X2_stream_vdc120);
    n4_s = not(and(n3_s, circshift(n1_s, 0)));%1
    n5_s = and(n4_s, X2_stream_vdc6);
    n6_s = not(and(n5_s, circshift(n1_s, 0)));%2
    y_s = and(n6_s, circshift(input,0));%2
    sin_sobol = sum(y_s)/N;


    input = X2_stream_lfsr;
    n1_lf = and(input, circshift(input,1));
    n2_lf = not(and(n1_lf, X2_stream_vdc5040));
    n3_lf = and(n2_lf, X2_stream_vdc120);
    n4_lf = not(and(n3_lf, circshift(n1_lf, 0)));%1
    n5_lf = and(n4_lf, X2_stream_vdc6);
    n6_lf = not(and(n5_lf, circshift(n1_lf, 0)));%2
    y_lf = and(n6_lf, circshift(input,0));%2
    sin_lfsr = sum(y_lf)/N;


end
end