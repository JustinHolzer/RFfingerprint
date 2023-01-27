function x_out = randRange(x, VecLen)
if(~exist('VecLen'))
    VecLen=1;
end
xmin = x(1); 
xmax = x(2);
x_out = (xmax-xmin).*rand(VecLen,1) + xmin;