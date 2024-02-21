function cos_lfsr = Cosine_lfsr_lfsr_Chu(X, N)
format long
if X > 1 || X < 0
    fprintf("Error");
else
%clear

%N = 1024;


if N == 1024
    [~,lfval] = LFSR_TrigonoSC([true false true true true false false false false false],X,N); %N=1024
    [~,lfval2] = LFSR_TrigonoSC2([false true false true false false false true false true],N/2,N); %N=1024
elseif N == 512
    [~,lfval] = LFSR_TrigonoSC([true false true false true false false false false],X,N); %N=512
    [~,lfval2] = LFSR_TrigonoSC2([true false true true false true false false true],N/2,N); %N=512
elseif N == 256
    [~,lfval] = LFSR_TrigonoSC([true false true false true false false false],N/2,N); %N=256
    [~,lfval2] = LFSR_TrigonoSC2([false false false false true true false false],N/2,N); %N=256
end

lfval = lfval/N;
lfval2 = lfval2/N;



%sin_vdc = zeros(1,N);


X2_stream_lfsr = zeros(1, N);

X2_stream_lfsr40320 = zeros(1, N);%1/56
X2_stream_lfsr720 = zeros(1, N);%1/30
X2_stream_lfsr24 = zeros(1, N);
X2_stream_lfsr2 = zeros(1, N);
%for i = 1:N
    for k = 1:N
        
        if X > lfval(k)
            X2_stream_lfsr(k) = 1;
        end
        if 1/40320 > lfval(k)
            X2_stream_lfsr40320(k) = 1;
        end
        if 1/720 > lfval(k)
            X2_stream_lfsr720(k) = 1;
        end
        if 1/24 > lfval(k)
            X2_stream_lfsr24(k) = 1;
        end
        if 1/2 > lfval(k)
            X2_stream_lfsr2(k) = 1;
        end
    end
%end

   
    

    n1_lf = and(X2_stream_lfsr, circshift(X2_stream_lfsr,1));
    n2_lf = not(and(n1_lf, circshift(X2_stream_lfsr40320,1)));
    n3_lf = and(n2_lf, circshift(X2_stream_lfsr720,1));
    n4_lf = not(and(n3_lf, circshift(n1_lf,1)));
    n5_lf = and(n4_lf, circshift(X2_stream_lfsr24,1));
    n6_lf = not(and(n5_lf, circshift(n1_lf,2)));
    n7_lf = and(n6_lf, circshift(X2_stream_lfsr2,1));
    y_lf = not(and(n7_lf, circshift(n1_lf,3)));
    cos_lfsr = sum(y_lf)/N;

end
end
