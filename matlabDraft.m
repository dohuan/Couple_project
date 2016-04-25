h_w_fitn4s2 = resultMatrix_h_w(end,:);
for i=1:7
	h_w_fit{i}=h_w_fitn4s2(1:i);
	h_w_fitn4s2(1:i)=[];
end

span = 3:30;
rmseTrain = zeros(size(span,1),1);
rmseTest = zeros(size(span,1),1);
count = 1;
for i = span
	setSplit = i;
	X_train = X(1:setSplit-1,:);
	X_test = X(setSplit:end,:);
	y_train = y(1:setSplit-1,:);
	y_test = y(setSplit:end,:);
	[FitObj,y_train_eval,y_guess_test]=lassoRegression(X_train,y_train,X_test,y_test);
	rmseTrain(count,1) = rmseCal(y_train_eval.y_guess_train,y_train);
	rmseTest(count,1) = rmseCal(y_guess_test,y_test);
	count=count+1;
end