function R2 = getR2(y,y_)
    SSres = sum((y-y_).^2);
    SStot = sum((y-mean(y)).^2);
    R2=1-SSres/SStot;
end