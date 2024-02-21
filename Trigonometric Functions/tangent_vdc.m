%Calculating the tan(x) using SC calculations.
%The input range is x is in [0,1] interval (in Radian)

function [tan_vdc, tan_sobol, tan_lfsr] = tangent_vdc(X, N)
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
X6_stream_sobol = zeros(1, N);
X2_stream_lfsr = zeros(1, N);

X2_stream_vdc24 = zeros(1, N);%1/42
X2_stream_vdc51 = zeros(1, N);%1/20
X2_stream_vdc170 = zeros(1, N);%1/6

X2_stream_vdc18 = zeros(1, N);%1/56
X3_stream_vdc34 = zeros(1, N);%1/30
X4_stream_vdc85 = zeros(1, N);%1/12
X5_stream_vdc512 = zeros(1, N);%1/2
%for i = 1:N
    for k = 1:N
        %if i/N > sobol(k,9)
        if X > vd(k,3)
            X2_stream_vdc(k) = 1;
        end
        if X > sobol(k,1)
            X2_stream_sobol(k) = 1;
        end
        if X > sobol(k,6)
            X6_stream_sobol(k) = 1;
        end
        if X > lfval(k)
            X2_stream_lfsr(k) = 1;
        end
        if 1/56 > vd(k,4)
            X2_stream_vdc18(k) = 1;
        end
        if 1/30 > vd(k,9)
            X3_stream_vdc34(k) = 1;
        end
        if 1/12 > vd(k,4)
            X4_stream_vdc85(k) = 1;
        end
        if 1/2 > vd(k,8)
            X5_stream_vdc512(k) = 1;
        end
        if 1/42 > vd(k,6)
            X2_stream_vdc24(k) = 1;
        end
        if 1/20 > vd(k,6)
            X2_stream_vdc51(k) = 1;
        end
        if 1/6 > vd(k,6)
            X2_stream_vdc170(k) = 1;
        end
    end
%end

%-------------------------Sine-------------------------------   
    n1_v = and(X6_stream_sobol, circshift(X6_stream_sobol,2));%best with 15 Delays
    n2_v = not(and(n1_v, X2_stream_vdc24));
    n3_v = not(and3(X2_stream_vdc51, n2_v, n1_v));
    n4_v = not(and3(X2_stream_vdc170, n3_v, n1_v));
    y_v = and(n4_v, X6_stream_sobol);
    sin_vdc = sum(y_v)/N;

    n1_s = and(X2_stream_sobol, circshift(X2_stream_sobol,2));%best with 15 Delays
    n2_s = not(and(n1_s, X2_stream_vdc24));
    n3_s = not(and3(X2_stream_vdc51, n2_s, n1_s));
    n4_s = not(and3(X2_stream_vdc170, n3_s, n1_s));
    y_s = and(n4_s, X2_stream_sobol);
    sin_sobol = sum(y_s)/N;

    n1_lf = and(X2_stream_lfsr, circshift(X2_stream_lfsr,2));%best with 15 Delays
    n2_lf = not(and(n1_lf, X2_stream_vdc24));
    n3_lf = not(and3(X2_stream_vdc51, n2_lf, n1_lf));
    n4_lf = not(and3(X2_stream_vdc170, n3_lf, n1_lf));
    y_lf = and(n4_lf, X2_stream_lfsr);
    sin_lfsr = sum(y_lf)/N;


%-------------------------Cosine-----------------------------
    n1_v_cos = and(X2_stream_vdc, circshift(X2_stream_vdc,2));%best with 15 Delays
    n2_v_cos = not(and(n1_v_cos, X2_stream_vdc18));
    n3_v_cos = not(and3(X3_stream_vdc34, n2_v_cos, n1_v_cos));
    n4_v_cos = not(and3(X4_stream_vdc85, n3_v_cos, n1_v_cos));
    y_v_cos = not(and3(X5_stream_vdc512, n4_v_cos, n1_v_cos));
    cos_vdc = sum(y_v_cos)/N;

    n1_s_cos = and(X2_stream_sobol, circshift(X2_stream_sobol,2));%best with 15 Delays
    n2_s_cos = not(and(n1_s_cos, X2_stream_vdc18));
    n3_s_cos = not(and3(X3_stream_vdc34, n2_s_cos, n1_s_cos));
    n4_s_cos = not(and3(X4_stream_vdc85, n3_s_cos, n1_s_cos));
    y_s_cos = not(and3(X5_stream_vdc512, n4_s_cos, n1_s_cos));
    cos_sobol = sum(y_s_cos)/N;

    n1_lf_cos = and(X2_stream_lfsr, circshift(X2_stream_lfsr,2));%best with 15 Delays
    n2_lf_cos = not(and(n1_lf_cos, X2_stream_vdc18));
    n3_lf_cos = not(and3(X3_stream_vdc34, n2_lf_cos, n1_lf_cos));
    n4_lf_cos = not(and3(X4_stream_vdc85, n3_lf_cos, n1_lf_cos));
    y_lf_cos = not(and3(X5_stream_vdc512, n4_lf_cos, n1_lf_cos));
    cos_lfsr = sum(y_lf_cos)/N;

 %--------------------------VDC-------------------------------------------------------------------
    if sin_vdc < cos_vdc
        [~,y_v_tan] = CORLD_DIV(N*sin_vdc, N*cos_vdc, N);
        tan_vdc = sum(y_v_tan)/N;
    else
        [sinus_vdc,~,~] = Sine_vdc(X/2,N);
        [cosinus_vdc,~,~] = Cosine_vdc(X/2,N);
        [~,y_v_tan] = CORLD_DIV(N*(sinus_vdc*cosinus_vdc), N*cos_vdc, N);
        %[~,y_v_tan] = CORLD_DIV(sum(and(y_v(floor(X/2)), y_v_cos(floor(X/2)))), N*cos_vdc, N);
        tan_vdc = 2*sum(y_v_tan)/N;
    end
%--------------------------Sobol-------------------------------------------------------------------
    if sin_sobol < cos_sobol
        [~,y_s_tan] = CORLD_DIV(N*sin_sobol, N*cos_sobol, N);
        tan_sobol = sum(y_s_tan)/N;
    else
        [~,sinus_sobol,~] = Sine_vdc(X/2,N);
        [~,cosinus_sobol,~] = Cosine_vdc(X/2,N);
        [~,y_s_tan] = CORLD_DIV(N*(sinus_sobol*cosinus_sobol), N*cos_sobol, N);
        %[~,y_s_tan] = CORLD_DIV(sum(and(y_s(floor(X/2)), y_s_cos(floor(X/2)))), N*cos_sobol, N);
        tan_sobol = 2*sum(y_s_tan)/N;
    end
%--------------------------LFSR-------------------------------------------------------------------
    if sin_lfsr < cos_lfsr
        [~,y_lf_tan] = CORLD_DIV(N*sin_lfsr, N*cos_lfsr, N);
        tan_lfsr = sum(y_lf_tan)/N;
    else
        [~,~,sinus_lfsr] = Sine_vdc(X/2,N);
        [~,~,cosinus_lfsr] = Cosine_vdc(X/2,N);
        [~,y_lf_tan] = CORLD_DIV(N*(sinus_lfsr*cosinus_lfsr), N*cos_sobol, N);
        %[~,y_lf_tan] = CORLD_DIV(sum(and(y_lf(floor(X/2)), y_lf_cos(floor(X/2)))), N*cos_lfsr, N);
        tan_lfsr = 2*sum(y_lf_tan)/N;
    end

end
end
