%Calculating the cos(x) using SC calculations.
%The input range is x is in [0,1] interval (in Radian)

function [cos_vdc, cos_sobol, cos_lfsr] = Cosine_vdc(X, N)
format long
if X > 1 || X < 0
    fprintf("Error");
else
%clear

%N = 1024;
sobol = net(sobolset(256), N);

[~,lfval] = LFSR_TrigonoSC([true false true true true false false false false false],X,N); %N=1024
%[~,lfval] = LFSR_TrigonoSC([true false true false true false false false false],X,N); %N=512
%[~,lfval] = LFSR_TrigonoSC([true false true false true false false false],X,N); %N=256
%[~,lfval] = LFSR_TrigonoSC([true false false true true false false],X,N); %N=128
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

X2_stream_vdc18 = zeros(1, N);
X3_stream_vdc34 = zeros(1, N);
X4_stream_vdc85 = zeros(1, N);
X5_stream_vdc512 = zeros(1, N);
%for i = 1:N
    for k = 1:N
        %if i/N > sobol(k,9)
        if X > vd(k,3)
            X2_stream_vdc(k) = 1;
        end
        if X > sobol(k,4)
            X2_stream_sobol(k) = 1;
        end
        if X > lfval(k)
            X2_stream_lfsr(k) = 1;
        end
        if 1/56 > vd(k,3)
            X2_stream_vdc18(k) = 1;
        end
        if 1/30 > vd(k,2)
            X3_stream_vdc34(k) = 1;
        end
        if 1/12 > vd(k,4)
            X4_stream_vdc85(k) = 1;
        end
        if 1/2 > vd(k,7)
            X5_stream_vdc512(k) = 1;
        end
    end
%end

   
    n1_v = and(X2_stream_vdc, circshift(X2_stream_vdc,2));%best with 15 Delays
    n2_v = not(and(n1_v, X2_stream_vdc18));
    n3_v = not(and3(X3_stream_vdc34, n2_v, n1_v));
    n4_v = not(and3(X4_stream_vdc85, n3_v, n1_v));
    y_v = not(and3(X5_stream_vdc512, n4_v, n1_v));
    cos_vdc = sum(y_v)/N;

    n1_s = and(X2_stream_sobol, circshift(X2_stream_sobol,2));%best with 15 Delays
    n2_s = not(and(n1_s, X2_stream_vdc18));
    n3_s = not(and3(X3_stream_vdc34, n2_s, n1_s));
    n4_s = not(and3(X4_stream_vdc85, n3_s, n1_s));
    y_s = not(and3(X5_stream_vdc512, n4_s, n1_s));
    cos_sobol = sum(y_s)/N;

    n1_lf = and(X2_stream_lfsr, circshift(X2_stream_lfsr,2));%best with 15 Delays
    n2_lf = not(and(n1_lf, X2_stream_vdc18));
    n3_lf = not(and3(X3_stream_vdc34, n2_lf, n1_lf));
    n4_lf = not(and3(X4_stream_vdc85, n3_lf, n1_lf));
    y_lf = not(and3(X5_stream_vdc512, n4_lf, n1_lf));
    cos_lfsr = sum(y_lf)/N;

end
end
