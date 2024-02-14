
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

%[~,lfval] = LFSR3([true false true true true false false false false false],N/2,N);
%LF3(LF3 == -1) = 0;
%lfval = lfval/N;
lfval = rand(1,N);

sq_root_sob = zeros(1,N);
abs_sq_root_sob = zeros(1,N);
sq_root_vd = zeros(1,N);
abs_sq_root_vd = zeros(1,N);
sq_root_lf = zeros(1,N);
abs_sq_root_lf = zeros(1,N);

EE = zeros(1,N);

X2_stream_sobol = zeros(N, N);
X2_stream_vdc = zeros(N, N);
X2_stream_lfsr = zeros(N, N);

for i = 1:N
    for k = 1:N
        if i/N > sobol(k,1)
            X2_stream_sobol(i,k) = 1;
        end
        if i/N > vd(k,2)
            X2_stream_vdc(i,k) = 1;
        end
        if i/N > lfval(k)
            X2_stream_lfsr(i,k) = 1;
        end
    end
end


for i = 1:N
    EE(i) = (i/N)^2;
    sq_root_sob(i) = sum(and(X2_stream_sobol(i,:), circshift(X2_stream_sobol(i,:),1)))/N;
    abs_sq_root_sob(i) = abs(EE(i) - sq_root_sob(i));
    sq_root_vd(i) = sum(and(X2_stream_vdc(i,:), circshift(X2_stream_vdc(i,:),1)))/N;
    abs_sq_root_vd(i) = abs(EE(i) - sq_root_vd(i));
    sq_root_lf(i) = sum(and(X2_stream_lfsr(i,:), circshift(X2_stream_lfsr(i,:),1)))/N;
    abs_sq_root_lf(i) = abs(EE(i) - sq_root_lf(i));
end
MAE_sobol = mean(abs_sq_root_sob)
%figure
%plot(abs_sq_root_sob)
%title('Sobol')
MAE_vdc = mean(abs_sq_root_vd)
%figure
%plot(abs_sq_root_vd)
%title('VDC')
MAE_lfsr = mean(abs_sq_root_lf)