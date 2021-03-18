function res = myrand(q)
    res = sum(rand >= cumsum(q));
end

