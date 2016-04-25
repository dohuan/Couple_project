function [ best_fit ] = MakeBestFitStruct( resultMatrix, max_value,max_RowIndex ,max_ColumnIndex)
% MakeBestFitStruct makes a struct for a best fit 
% it includes the best fit value, model number, split position, start index
% ,end index for the input, and the width of the split.

    best_fit = struct('value',max_value);
    best_fit.modelNumber =max_RowIndex;
    best_fit.splitNumber = resultMatrix(1,max_ColumnIndex(max_RowIndex));
    best_fit.startPos = resultMatrix(2,max_ColumnIndex(max_RowIndex));
    best_fit.endPos = resultMatrix(3,max_ColumnIndex(max_RowIndex));
    best_fit.width = best_fit.endPos -best_fit.startPos+1;
    
    
end

