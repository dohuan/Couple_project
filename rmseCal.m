function rmse = rmseCal(predict,target)
    [m,~] = size(predict);
    err_tmp = predict - target;
    err = sum(err_tmp.^2,2);
    rmse = sqrt(sum(err)/m);
end