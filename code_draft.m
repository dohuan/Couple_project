subplot(3,1,1)
plot(y)
subplot(3,1,2)
hold on
for i=size(X,2)-5:size(X,2)
    plot(X(:,i),'--');
end
hold off
subplot(3,1,3)
plot(X_mu)

for i=1:size(X,2)
	cov_matrix(:,:,i) = cov([zscore(y) zscore(log(X(:,i)))]);
end
 % mode 3 and 4 have highest correlation with y: 0.248 and 0.2786
 
 
 
 % -- use LASSO to find out useful features
 [ix,ic] = find(abs(FitObj.B_optimal)>.05);
 fweight = FitObj.B_optimal(ix);
 for i=1:length(ix)
	figure(i)
	hold on
	plot(y);
	plot(X(:,i));
	hold off
	cov([y X(:,i)])
 end
 
 y_ = X(:,ix)*fweight;
 
 % --- Plots for ppt
 y_ = X(:,ix)*FitObj.B_optimal(FitObj.);