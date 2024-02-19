function cos_lfsr = Cosine_lfsr_lfsr_Parhi(X, N)
format long
if X > 1 || X < 0
    fprintf("Error");
else
%clear

%N = 1024;


[~,lfval] = LFSR_TrigonoSC([true false true true true false false false false false],N/2,N); %N=1024
%[~,lfval] = LFSR_TrigonoSC([true false true false true false false false false],N/2,N); %N=512
%[~,lfval] = LFSR_TrigonoSC([true false true false true false false false],N/2,N); %N=256
%[~,lfval] = LFSR_TrigonoSC([true false false true true false false],N/2,N); %N=128
lfval = lfval/N;



%sin_vdc = zeros(1,N);


X2_stream_lfsr = zeros(1, N);

X2_stream_lfsr56 = zeros(1, N);%1/56
X2_stream_lfsr30 = zeros(1, N);%1/30
X2_stream_lfsr12 = zeros(1, N);
X2_stream_lfsr2 = zeros(1, N);
%for i = 1:N
    for k = 1:N
        
        if X > lfval(k)
            X2_stream_lfsr(k) = 1;
        end
        if 1/56 > lfval(k)
            X2_stream_lfsr56(k) = 1;
        end
        if 1/30 > lfval(k)
            X2_stream_lfsr30(k) = 1;
        end
        if 1/12 > lfval(k)
            X2_stream_lfsr12(k) = 1;
        end
        if 1/2 > lfval(k)
            X2_stream_lfsr2(k) = 1;
        end
    end
%end

   
    

    n1_lf = and(X2_stream_lfsr, circshift(X2_stream_lfsr,4));
    n2_lf = not(and(n1_lf, circshift(X2_stream_lfsr56,1)));
    n3_lf = not(and3(circshift(X2_stream_lfsr30,1), n2_lf, circshift(n1_lf,1)));
    n4_lf = not(and3(circshift(X2_stream_lfsr12,1), n3_lf, circshift(n1_lf,2)));
    y_lf = not(and3(circshift(X2_stream_lfsr2,1), n4_lf, circshift(n1_lf,3)));
    cos_lfsr = sum(y_lf)/N;

end
end