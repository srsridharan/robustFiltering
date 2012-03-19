function quadcell = fitForFn(fntofit, xhat, numberOfQuadratics, numOfPtsAroundEachQuadratic, rangeAroundXhat )

% global GLOBAL_optimsetting;

	% opt1= optimset(GLOBAL_optimsetting,'Display','off');	
  xhat = xhat(2);
	opt1= optimset('Display','off');	
	
	quadcell = {};
%{
	numberOfQuadratics = 7;
	numOfPtsAroundEachQuadratic = 20;
	rangeAroundXhat = 3;
	%}

%  keyboard;
	xsamplePts = linspace(xhat - rangeAroundXhat, xhat + rangeAroundXhat, 30)';
	fnvalatxsamplepts = fntofit(xsamplePts);
% 	A = [ones(length(xsamplePts),1), xsamplePts, xsamplePts.^2];
% 	fnewsamplepts= - 2*yval'*A;
% 	Hsamplepts = 2*A'*A;
		
	centerPts = linspace(xhat-rangeAroundXhat, xhat + rangeAroundXhat, numberOfQuadratics + 1);
%     keyboard;
	for(k = 1:length(centerPts)-1)
		xval = linspace(centerPts(k),centerPts(k+1),numOfPtsAroundEachQuadratic)';
		yval = fntofit(xval);
		A = [xval.^2,xval, ones(length(xval),1)];
		f = -2*A'*yval;
		fnew= f';
		H = 2*A'*A;
		% quadprog minimizes (1/2)*x'*H*x + f'*x
% 		pval = fmincon(@(x) ((1/2)*x'*H*x + fnew*x), [150;80;160],[],[],[],[],[1;-250;-250],[],@(x)fnnonlincon(x,xsamplePts,fnvalatxsamplepts));
        

% for + - x^3/40
         pval = fmincon(@(x) ((1/2)*x'*H*x + fnew*x), [3/40;mean(diff(yval))/mean(diff(xval));mean(yval)],[],[],[],[],[.01;-250;-250],[],@(x)fnnonlincon(x,xsamplePts,fnvalatxsamplepts),opt1);

        
        
        
        
     % for x^6/40^3  
%      ver 1 : second derivative: 6*xhat^3/(40^2
%        pval = fmincon(@(x) ((1/2)*x'*H*x + fnew*x), [(30/40^2)*(mean(xval))^4;mean(diff(yval)./diff(xval));fntofit(mean(xval))],[],[],[],[],[.01;-250;-250],[],@(x)fnnonlincon(x,xsamplePts,fnvalatxsamplepts));
        %   =                 fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
        
        %{
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
		pval = quadprog(H,f,[],[],[],[],[-200;-200;0],[]);
		Q = [0 0 0;...
				 0 pval(3) pval(2)/2;...
				 0 pval(2)/2 pval(1)];
%}
		quadcell{k} = [0 0 0;...
									0 pval(1) pval(2)/2;...
									0 pval(2)/2 pval(3)];


		%{
		figure;
        plot(xval,fntofit(xval),'b--');
        hold on ;
        yvalmin = polyval(quadcell{k},xval);
        plot(xval, yvalmin,'r-.');
	%}	



	end


end


function [c,ceq] = fnnonlincon(x,xsamplePts,fnvalatxsamplepts)
% enhancements: return the gradients analytically
		ceq = 0;
%         keyboard;
		c = max(fnvalatxsamplepts - polyval(x,xsamplePts));



end
