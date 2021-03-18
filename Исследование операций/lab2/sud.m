function res = sud(A,B)

if (A(2)==B(1) && B(2)==A(1))
    res = 0;
elseif  (A(2)~=B(1) && B(2)~=A(1))
    res = 0;
elseif (A(2)==B(1) && B(2)~=A(1))
    res = sum(A);
else
    res = -sum(B);
end
end