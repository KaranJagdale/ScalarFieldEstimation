function [bis] = perpbisector(line1)
% this function computes the perpendicular bisector to 'line1' - the perpendicular bisector is formed by the line connecting points (xnew1,ynew1) and (xnew2,ynew2)

	TOL = 0.0000001;

	% compute mid-point of 'line1' - this is a point on the perpendicular bisector
	xnew1 = (line1.x1 + line1.x2)/2;
	ynew1 = (line1.y1 + line1.y2)/2;
	
	% slope of the perpendicular bisector is close to zero
	if(((line1.x2 - line1.x1)<TOL) && ((line1.x2 - line1.x1)>-TOL))
		m = 0;
		% choose a point with x-coordinate shifted by 0.5.. this is for finding another point on the perpendicular bisector
		xnew2 = xnew1 + 0.5;
		% since slope~0, the y-coordinates do not change.. 
		ynew2 = ynew1;
	else
		% slope is zero for line1 and consequently the perp bisector has slope infinity
		if(((line1.y2 - line1.y1)<TOL) && ((line1.y2 - line1.y1)>-TOL))
			ynew2 = ynew1 + 0.5;
			xnew2 = xnew1;
		else
		% the case where 'line1' has neither slope zero nor slope infinity
			% compute the slope of 'line1'
			m1 = (line1.y2 - line1.y1)/(line1.x2 - line1.x1);
			% slope of the perpendicular bisector
			m = -1/m1;
			% if the slope of perp bisector is such that it is more aligned with the x axis
			if((m<=1.0) && (m>=-1.0))
				xnew2 = xnew1 + 0.5;
				ynew2 = (xnew1*(line1.x1-line1.x2) + ynew1*(line1.y1-line1.y2) - xnew2*(line1.x1-line1.x2))/(line1.y1-line1.y2);
			else
			% if the slope of perp bisector is such that it is more aligned with the y axis
				ynew2 = ynew1 + 0.5;
				xnew2 = ((ynew2-ynew1) + m*xnew1)/m;
			end
		end
	end
	
	% store the perpendicular bisector constructed from the points (xnew1,ynew1) and (xnew2,ynew2) in the structure 'bis'
	bis = struct('x1',xnew1,'y1',ynew1,'x2',xnew2,'y2',ynew2);

end