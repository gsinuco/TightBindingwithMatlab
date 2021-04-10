function [T] = translation()
    global L

    T = 0.0;

    for i=1:L-1;
        T(i,i+1) = 1.0;
        T(i+1,i) = 1.0;
    end

    T(L,1) = 1.0;
    T(1,L) = 1.0;