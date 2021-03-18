function res = treug1(x,v,tau0,tau1,tau2)
Qm = v * tau2;
lamb = Qm / tau1 + v;
if ((x <= tau0) || (x >= tau0 + tau1 + tau2))
    res = 0;
else if ((x <= tau0 + tau1) && (x >= tau0))        
        res = (lamb - v).*(x - tau0);
    else if ((x <= tau0 + tau1 + tau2) && (x >= tau1))
            res = Qm - v.*(x - (tau1 + tau0));
        end
    end
end
end