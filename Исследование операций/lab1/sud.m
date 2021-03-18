function [AV, BV] = sud(A,B,AV,BV)
if A(2)==B(1) & B(2)==A(1)
    return
elseif    A(2)~=B(1) & B(2)~=A(1)
    return
elseif A(2)==B(1) & B(2)~=A(1)
    AV=AV+sum(A);
else
    BV=BV+sum(B);
end    
%[AV BV]    
end