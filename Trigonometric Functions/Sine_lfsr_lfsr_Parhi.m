function sin_lfsr = Sine_lfsr_lfsr_Parhi(X, N)
format long
if X > 1 || X < 0
    fprintf("Error");
else
%clear

%N = 1024;


[~,lfval] = LFSR_TrigonoSC([true false true true true false false false false false],N/2,N); %N=1024
%[~,lfval] = LFSR_TrigonoSC([true false true false true false false false false],N/2,N); %N=512
lfval = lfval/N;



%sin_vdc = zeros(1,N);


X2_stream_lfsr = zeros(1, N);
X2_stream_lfsr24 = zeros(1, N);
X2_stream_lfsr51 = zeros(1, N);
X2_stream_lfsr170 = zeros(1, N);


%for i = 1:N
    for k = 1:N
        %if i/N > sobol(k,9)
        
        if X > lfval(k)
            X2_stream_lfsr(k) = 1;
        end
        if 1/42 > lfval(k)
            X2_stream_lfsr24(k) = 1;
        end
        if 1/20 > lfval(k)
            X2_stream_lfsr51(k) = 1;
        end
        if 1/6 > lfval(k)
            X2_stream_lfsr170(k) = 1;
        end
        
    end
%end

    

    n1_lf = and(X2_stream_lfsr, circshift(X2_stream_lfsr,3));
    n2_lf = not(and(n1_lf, circshift(X2_stream_lfsr24,1)));
    n3_lf = not(and3(circshift(X2_stream_lfsr51,1), n2_lf, circshift(n1_lf,1)));
    n4_lf = not(and3(circshift(X2_stream_lfsr170,1), n3_lf, circshift(n1_lf,2)));
    y_lf = and(n4_lf, circshift(X2_stream_lfsr,6));
    sin_lfsr = sum(y_lf)/N;

end
end