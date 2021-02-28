function [ynew, tnew] = eulermethod(told, yold, h)
    ynew = yold + h*yold;
    tnew = told + h;
end