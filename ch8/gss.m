function [fmax,argmax] = gss(xmin,xmax,fn)
% Golden section search to find the max of a single peaked function
p = 0.618033988749895;
A = xmin; D = xmax;
B = p*A + (1-p)*D;
C = (1-p)*A + p*D;
f_B = fn(B); f_C = fn(C);

epsilon = 1e-3;
dist = 10;
while dist > epsilon*max(1,abs(B)+abs(C))
    if f_B > f_C
        D = C; C = B;
        f_C = f_B;
        B = p*C+(1-p)*A;
        f_B = fn(B);
    else
        A = B; B = C;
        f_B = f_C;
        C = p*B+(1-p)*D;
        f_C = fn(C);
    end
    dist = abs(D-A);
end
argmax = B;
fmax = f_B;
end