function [xAxis,yAxis,sigCurveParams] = fitSigmoidToData(xAxis,yAxis, xStart,xStop, initValue)

    % We define the parameters for the fit
    fo = fitoptions('Method', 'NonlinearLeastSquares', ...
    'StartPoint', initValue, ...
    'Lower', [-inf, -inf, 0, 0, -inf, -inf], ...
    'MaxFunEvals', 12000000, 'MaxIter', 800000000, ...
    'TolFun', 1e-12);
    ft = fittype('a + b./(c + d.*exp(-f.*(x - g)))', ...
    'options', fo);

   % we exclude data from the fit 
     yAxis = yAxis(find(xAxis < xStop));
     xAxis = xAxis(find(xAxis < xStop));
     yAxis = yAxis(find(xAxis > xStart));
     xAxis = xAxis(find(xAxis > xStart));

   % We calculate the fit
    sigCurveParams = fit(xAxis, yAxis, ft);

end