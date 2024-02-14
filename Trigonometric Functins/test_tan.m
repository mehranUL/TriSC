clear
X = 0:0.01:1;
N = 1024;
y_vdc = zeros(1,length(X));
y_sobol = zeros(1,length(X));
y_lfsr = zeros(1,length(X));

for i = 1:length(X)
    [y_vdc(i),y_sobol(i),y_lfsr(i)] = tan_vdc(X(i), N);
end

figure
plot(y_vdc,'-*')
hold on
plot(y_sobol,'-o')
plot(y_lfsr,'-+')
plot(tan(X),'-^')