function [  ] = PlotFits( result , wife, husband, time)
%PlotFit plot the given fit
    
    
    wife = wife - mean(wife);
    husband = husband - mean(husband(:,:));
    
    figure('name','basic')
    hold on
    plot(time(1:size(wife,1)),wife,time(1:size(wife,1)),husband)
    legend('wife','husband')

    startP=result.best_fit_w.startPos;
    endP=result.best_fit_w.endPos;
    wife1 = wife(startP:endP,:);
    husband1 = husband(startP:endP,:);
    
    figure('name','best fit for input:wife, model:123)')
    hold on
    plot(time(startP:endP),wife1,time(startP:endP),husband1)
    legend('wife','husband')
    
    
    startP2=result.best_fit_h.startPos;
    endP2=result.best_fit_h.endPos;
    wife2 = wife(startP2:endP2,:);
    husband2 = husband(startP2:endP2,:);
    
    figure('name','best fit for input:husband, model:123)')
    hold on
    plot(time(startP2:endP2),wife2,time(startP2:endP2),husband2)
    legend('wife','husband')
     

end

