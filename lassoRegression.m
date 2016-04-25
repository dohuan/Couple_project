function [FitObj,y_train_eval,y_test_eval]=lassoRegression(XTrain,yTrain,XTest,yTest,validate)
    %----------------------------------------------------------------------
    % -- Apply LASSO to predict output with given train and validate set --
    %                      ~~ Author: Huan N. Do ~~
    %----------------------------------------------------------------------
    % [FitObj,y_train_eval,y_test_eval]=
    %                   lassoRegression(XTrain,yTrain,XTest,yTest,validate)
    % validate: radio of data used for training over whole train set
    %                       (no input for DEactivate)
    % FitInfo: fitting info from LASSO 
    % mode = 0: use Test data to find optimal lambda
    % mode = 1: use validate set to find optimal lambda
    % y_train_eval.{y_guess_train,y_train}
    %----------------------------------------------------------------------
    if nargin <4
        error('myApp:argChk', 'Not enough inputs!')
    else if nargin <5
            mode = 0; %use Test data to find optimal lambda
        else
            mode = 1; % use validate set to find optimal lambda
            if validate<0||validate>1 
                error('myApp:argChk', 'Wrong train-validate ratio!')
            end
        end
    end
    
    x_test = XTest;
    y_test = yTest;
    testSize = size(x_test,1);
    if mode == 1
        splitIndex = round(size(XTrain,1)*validate);
        x_val = XTrain(splitIndex:end,:);
        y_val = yTrain(splitIndex:end,:);

        x_train = XTrain(1:splitIndex-1,:);
        y_train = yTrain(1:splitIndex-1,:);
    else 
        x_train = XTrain;
        y_train = yTrain;
    end 
    
    %trainSize = size(x_train,1);
    [B,FitInfo]=lasso(x_train,y_train);
    total=length(FitInfo.Intercept); %total number of lambda's
    if mode == 1
        valSize = size(x_val,1);
        y_guess = zeros(valSize,total);
        for i=1:valSize
            for j=1:total
                y_guess(i,j)=x_val(i,:)*B(:,j)+ FitInfo.Intercept(j);
            end
        end
        val_error = zeros(total,1);
        for i=1:total
            val_error(i)=mean((y_val-y_guess(:,i)).^2)/total;
        end
        minIndex=find(val_error==min(val_error));
        BOpt=B(:,minIndex); %column of B optimally results smallest MSE
        interceptOpt=FitInfo.Intercept(minIndex);
        s=sprintf('min index of lambda: %d\n',minIndex);
        disp(s);
        %----------------use optimal B to predict the TEST set-------------
        y_test_eval = x_test*BOpt+interceptOpt;
        y_error = y_test - y_test_eval;
        mse = sum(y_error.^2)/size(y_error,1);
        %----------------use optimal B to predict the TRAIN set------------
        y_train_eval.y_guess_train = x_train*BOpt+interceptOpt;
        y_train_eval.y_train = y_train;
        %--------------------------save Fit object-------------------------
        FitObj.Fit = FitInfo;
        FitObj.B = B;
        FitObj.minIndex = minIndex;
        FitObj.B_optimal = BOpt;
        FitObj.intercept_optimal = interceptOpt;
        FitObj.mse = mse;
        FitObj.splitIndex = splitIndex;
        FitObj.validate_error_curve = val_error;
        
    else 
        y_guess = zeros(testSize,total);
        for i=1:testSize
            for j=1:total
                y_guess(i,j)=x_test(i,:)*B(:,j)+ FitInfo.Intercept(j);
            end
        end
        val_error = zeros(total,1);
        for i=1:total
            val_error(i)=mean((y_test-y_guess(:,i)).^2)/total;
        end
        minIndex=find(val_error==min(val_error));
        BOpt=B(:,minIndex); %column of B optimally results smallest MSE
        interceptOpt=FitInfo.Intercept(minIndex);
        s=sprintf('min index of lambda: %d\n',minIndex);
        disp(s);
        %----------------use optimal B to predict the TEST set-------------
        y_test_eval = x_test*BOpt+interceptOpt;
        y_error = y_test - y_test_eval;
        mse = sum(y_error.^2)/size(y_error,1);
        %----------------use optimal B to predict the TRAIN set------------
        y_train_eval.y_guess_train = x_train*BOpt+interceptOpt;
        y_train_eval.y_train = y_train;
        %--------------------------save Fit object-------------------------
        FitObj.Fit = FitInfo;
        FitObj.B = B;
        FitObj.minIndex = minIndex;
        FitObj.B_optimal = BOpt;
        FitObj.intercept_optimal = interceptOpt;
        FitObj.mse = mse;
        FitObj.test_error_curve = val_error;
    end
    
end