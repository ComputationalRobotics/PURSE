 function [R,t] = solveP3P(y,Y)
    y     = normc([y;ones(1,3)]);
    [R,t] = p3p_nakano_bmvc2019(y,Y);
end