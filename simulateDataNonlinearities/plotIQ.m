function plotIQ(x, figNum,plotColor)
xi = real(x);
xq = imag(x);
if(~exist('figNum'))
    figure
else
    figure(figNum)
end
plot(xi,xq,'x')
xlabel('in-phase')
ylabel('quadrature')