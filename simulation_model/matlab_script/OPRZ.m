function J = OPRZ(vhI,vhP);
P2 = rand;
J = 0;
B = 0;
x = size(vhP);
while  J < x(2)
    if (P2>B)
        J = J + 1;
        B = B + vhP(vhI,J);     
    else
        break;
    end
end

