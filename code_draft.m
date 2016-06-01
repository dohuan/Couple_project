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
 
% --- Plot weight with name
loadFile{1}= 'couple_both_30_warm';
loadFile{2}= 'couple_both_30_domi';
loadFile{3}= 'couple_husband_warm';
loadFile{4}= 'couple_wife_warm';
loadFile{5}= 'couple_husband_domi';
loadFile{6}= 'couple_wife_domi';
savePath = 'C:\Users\dohuan\OneDrive\Graduate Research(DB)\CoupleProject\presentationForChris\';

for i=1:length(loadFile)
	load(loadFile{i});
	if (i==1)
		f_name = couple(1).f_name;
	end
	h = figure(1);
	bar(FitObj.B_optimal);
	ix = find(abs(FitObj.B_optimal)>0.5);
	B_select = FitObj.B_optimal(ix);
	for j=1:length(ix)
		text(ix(j)+2, B_select(j), f_name{ix(j)},'BackgroundColor',[1 1 1]);
	end
	
	axis tight
	box on
	ylabel('weights')
	xlabel('features')
	saveas(h,[savePath loadFile{i} '_fname.jpg']);
	
end
% --- Compute R squared

loadFile{1}= 'couple_both_30_warm';
loadFile{2}= 'couple_both_30_domi';
loadFile{3}= 'couple_husband_warm';
loadFile{4}= 'couple_wife_warm';
loadFile{5}= 'couple_husband_domi';
loadFile{6}= 'couple_wife_domi';
for i=1:length(loadFile)
	load(loadFile{i});
	fprintf(['R2 ' loadFile{i} ':%f\n'],getR2(y,y_fit));
	
end


% --- Plot stat of eigenvalue
%eigIx = [2:6 23:27 44:48 65:69 86:90];

loadFeature = 'warm';
loadFile = ['couple_both_' loadFeature '_CV_H2.mat'];
load(loadFile);
close all
eigIx = [2:6];
X_ = [];
for i=1:kfold
	test_ix = index(i+(i-1)*(spf-1):i+i*(spf-1));
	for j=1:length(test_ix)
			X_ = [X_; couple(test_ix(j)).feature(1,eigIx)];
	end
	%X_ = [X_; X(test_ix,eigIx)];
end
%for i=1:length(couple)
%	X_ = [X_;couple(i).feature(1,eigIx)];
%end
figure(1)
suptitle(['Both ' loadFeature ' H2'])
for i=1:length(eigIx)
	subplot(1,5,i)
	hist(X_(:,i),10);
	xlabel(['eig-' num2str(i)]);
	ylim([0 20])
end


% --- check descending eigvalue
eigIx = [2 19 36 53 70];

