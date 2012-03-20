function quadcell = fitForFnSq(fntofit, xhat, numberOfQuadratics, numOfPtsAroundEachQuadratic, rangeAroundXhat )
%% this function returns a min-plus expansion of the square of the function fntofit around xhat. It does this by using a 
% nonlinear constrained optimization on a grid of numOfPtsAroundEachQuadratic  around the center point of each
% quadratic. The centerpoints of each quadratic are chosen by linearly choosing points in a window zone (rangeAroundXhat) 
% about  xhat. 

	opt1= optimset('Display','off');% it speeds up the run if messages arent being printed out. Comment this line during debugging		
	quadcell = {};
  xhat = xhat(2); % in the example used only the 2nd state influences the output directly (since y = ((x_2) ^3)/40

	xsamplePts = linspace(xhat - rangeAroundXhat, xhat + rangeAroundXhat, 30)';
	fnvalatxsamplepts = fntofit(xsamplePts);
	centerPts = linspace(xhat-rangeAroundXhat, xhat + rangeAroundXhat, numberOfQuadratics + 1);
	for(k = 1:length(centerPts)-1)
		xval = linspace(centerPts(k),centerPts(k+1),numOfPtsAroundEachQuadratic)';
		yval = fntofit(xval);
		A = [xval.^2,xval, ones(length(xval),1)];
		f = -2*A'*yval;
		fnew= f';
		H = 2*A'*A;
		pval = fmincon(@(x) ((1/2)*x'*H*x + fnew*x), [(30/40^2)*(mean(xval))^4;mean(diff(yval)./diff(xval));fntofit(mean(xval))],...
					[],[],[],[],[.01;-250;-250],[],@(x)fnnonlincon(x,xsamplePts,fnvalatxsamplepts),opt1);
		quadcell{k} = [0 0 0;...
									0 pval(1) pval(2)/2;...
									0 pval(2)/2 pval(3)];
	end % of for loop


end % of fn fitForFnSq


function [c,ceq] = fnnonlincon(x,xsamplePts,fnvalatxsamplepts)
% the nonlinear constraint is that the quadratic must always lie above the function being fitted in the zone of interst.
% TODO: enhancements: return the gradients analytically

	ceq = 0;
		c = max(fnvalatxsamplepts - polyval(x,xsamplePts));

end % of fn fnnonlincon
