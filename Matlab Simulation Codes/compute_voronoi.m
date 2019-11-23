%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLAB code for computing the voronoi partition of the j-th site in a bounded polygon %%
%% Author: Rihab Abdul Razak, IITB-Monash Research Academy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ yvert,xvert] = compute_voronoi(j, xborder, yborder, xsites, ysites)
% xborder - the x-coordinates of the boundary
% yborder - the y-coordinates of the boundary
% xsites - the x-coordinates of the sites
% ysites - the y-coordinates of the sites
% xvert - gives the x-coordinates of the voronoi cell of site j
% yvert - gives the y-coordinates of the voronoi cell of site j
% loop this function over j to get all the voronoi cells

	% Get the coordinates for my site (j-th site)
	myx = xsites(j);
	myy = ysites(j);

	n = length(xsites);
	if(length(ysites) ~= n)
		disp('The length of xsites and ysites must be same..');
		xvert = [];
		yvert = [];
		return;
	end

	% new variables for storing sorted x and y coordinates
	sortedxsites = xsites;
	sortedysites = ysites;

	% erase the j-th element from site list
	sortedxsites(j) = [];
	sortedysites(j) = [];

	% sort elements based on distance to site j
	for i=1:n-1
		distancetoagents(i) = (myx-sortedxsites(i))^2 + (myy-sortedysites(i))^2;
	end
	[sortedxsites, sortedysites] = sortsites(j, n, distancetoagents, sortedxsites, sortedysites);

	borderx = xborder;
	bordery = yborder;

	% curr and bis are structures which represent a line segment between (x1,y1) and (x2,y2)
	curr = struct('x1',[],'y1',[],'x2',[],'y2',[]);
	%bis = struct('x1',[],'y1',[],'x2',[],'y2',[]);

	% add my site as the starting point of the line 'curr'
	curr.x1 = myx;
	curr.y1 = myy;

	% loop over the sites other than j
	for i=1:n-1
	
		% get the i-th site coordinates
		xcur = sortedxsites(i);
		ycur = sortedysites(i);

		% add the i-th site as the second point defining the line 'curr'
		curr.x2 = xcur;
		curr.y2 = ycur;
		% compute the perpendicular bisector of line 'curr'
		[bis] = perpbisector(curr);
		% compute the intersection of the perpendicular bisector (line 'bis') with the polygon defined by borderx and bordery and store the resulting polygon (which lies on the same side of the line as that of my site) in the same variables borderx and bordery
		[borderx, bordery] = lineintersection(borderx, bordery, bis, myx, myy);
	
	end

	% the required voronoi partition is given by xvert and yvert
	xvert = borderx;
	yvert = bordery;

end

function [sx, sy] = sortsites(j, n, distancetoagents, sortedxsites, sortedysites)
% sort the sites based on distancetoagents

	sx = sortedxsites;
	sy = sortedysites;

	% sortedxsites and sortedysites are of length (n-1)
	for i=1:(n-2)
		for k=(i+1):(n-1)
			if(distancetoagents(i)>distancetoagents(k))
				temp = distancetoagents(i);
				distancetoagents(i) = distancetoagents(k);
				distancetoagents(k) = temp;

				temp = sx(i);
				sx(i) = sx(k);
				sx(k) = temp;

				temp = sy(i);
				sy(i) = sy(k);
				sy(k) = temp;
			end
		end
	end

end


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

% intersection of a line with polygon defined by its vertices
function [borderx, bordery] = lineintersection(borderx, bordery, line, myx, myy)
% this function returns a new polygon formed by the intersection of a line with a given polygon.. the new polygon is that which lies on that side of the line which contains the point (myx,myy)

	TOL = 0.0000001;
	condition = 0;
	% determine if the slope of the line is infinity - condition 1
	if(((line.x2 - line.x1)<TOL) && ((line.x2 - line.x1)>-TOL)) 
		condition = 1;
	end
	
	% case where the slope of line is infinity
	if(condition==1) 
		% check if myx is larger than the x-coordinate of the line.. if it is larger, then it is on the right side of the line (side = 1); otherwise it is on the left side of the line (side = -1).
		if((myx-line.x1)>0) 
			side = 1;
		else 
			side = -1;
		end
	% case where the slope of the line is not infinity
	else
	
		% compute the slope of the line
		m = (line.y2 - line.y1)/(line.x2 - line.x1);

		% determine which side of the line the point (myx,myy) is
		if(((myy-line.y1) - m*(myx-line.x1)) > 0) 
			side = 1;
		else 
			side = -1;
		end
	end

	% std::vector<double> borderxtmp(borderx);
	% std::vector<double> borderytmp(bordery);
	borderxtmp = borderx;
	borderytmp = bordery;

	% n = static_cast<int>(borderx.size());
	n = length(borderx);
	% determine the polygon vertices which lie on the same side as the point
	% std::vector<int> ind(n,0);
	ind = zeros(n,1);

	% to track if no changes have been done to the input polygon
	nochange = 1;

	% loop over the vertices of the polygon
	for i=1:n
	
		% the current vertex
		currx = borderx(i);
		curry = bordery(i);

		switch side	
			% if side = 1
			case 1 
				% check if the current vertex i is also on the same side of the line.. if it so set ind(i) = 1
				if(condition==1)  
					if((currx-line.x1)>=0) 
						ind(i) = 1;
					else 
						nochange = 0;
					end
				else 
					if(((curry-line.y1) - m*(currx-line.x1)) >= 0) 
						ind(i) = 1;
					else  
						nochange = 0;
					end
				end

			case -1 
				% check if the current vertex i is also on the same side of the line.. if it so set ind(i) = 1
				if(condition==1) 
				
					if((currx-line.x1)<=0) 
						ind(i) = 1;
					else 
						nochange = 0;
					end
				else 
					if(((curry-line.y1) - m*(currx-line.x1)) <= 0) 
						ind(i) = 1;
					else 
						nochange = 0;
					end
				end
		end
	end

	% if no changes to be done, return from the function
	if(nochange==1)
		return
	end

	%int i1,i2;
	%double ix, iy;
	%struct linedesc line2;

	% new line 'line2' to represent a segment of the boundary of the input polygon which needs to be cut
	line2 = struct('x1',[],'y1',[],'x2',[],'y2',[]);
	i1 = 0;
	i2 = 0;

	% if the first vertex is on the wrong side of the line, it needs to be removed and replaced
	if(ind(1)==0)
	
		% loop over vertices 1 to n-1
		for i=1:n-1
		
			if((ind(i)==0) && (ind(i+1)==1)) 
			
				% all vertices upto vertex i are on the wrong side - store i as i1
				i1 = i;

				% construct line2 joining the vertices i and i+1
				line2.x1 = borderx(i);
				line2.y1 = bordery(i);
				line2.x2 = borderx(i+1);
				line2.y2 = bordery(i+1);
				% compute intersection
				[ix,iy] = computeintersection(line, line2);
				
				% replace and update the boundary vertex i with the intersection (ix,iy)
				borderxtmp(i) = ix;
				borderytmp(i) = iy;
			end
			if((ind(i)==1) && (ind(i+1)==0))
			
				% all vertices starting from vertex i+1 are on the wrong side - store i+1 as i2
				i2 = i+1;

				% construct line2 joining the vertices i and i+1
				line2.x1 = borderx(i);
				line2.y1 = bordery(i);
				line2.x2 = borderx(i+1);
				line2.y2 = bordery(i+1);
				% compute intersection
				[ix,iy] = computeintersection(line, line2);
				
				% replace and update the boundary vertex i+1 with the intersection (ix,iy)
				borderxtmp(i+1) = ix;
				borderytmp(i+1) = iy;
			end
		end
		% if the last vertex is on the same side of the line
		if(ind(n)==1) 
		
			% set i2=1.. the set of vertices on the wrong side of the line is from i2 to i1..
			i2 = 1;
	
			% construct line2 joining the vertices n and 1
			line2.x1 = borderx(n);
			line2.y1 = bordery(n);
			line2.x2 = borderx(1);
			line2.y2 = bordery(1);
			% compute intersection
			[ix,iy] = computeintersection(line, line2);
				
			% update the boundary vertices
			if(i1==i2)
			% only the first vertex is on the wrong side..
            %then i need to delete the vertex and replace 
            %it by two new vertices which are the intersections of
            %the line with segment of border joining the first vertex 
            %with the neighbouring vertices.. the first vertex is replaced by 
            %one new vertex already.. i need to add the other new vertex at the end..
				%borderxtmp.push_back(ix);
				%borderytmp.push_back(iy);
				borderxtmp(end+1) = ix;
				borderytmp(end+1) = iy;
			else
			% there are atleast two vertices on the wrong side..
				borderxtmp(1) = ix;
				borderytmp(1) = iy;
			end

			if((i2+1)<(i1)) 
			%borderxtmp.erase(borderxtmp.begin() + (i2 + 1), borderxtmp.begin() + (i1));
			%borderytmp.erase(borderytmp.begin() + (i2 + 1), borderytmp.begin() + (i1));
			borderxtmp(i2+1:i1-1) = [];
			borderytmp(i2+1:i1-1) = [];
			end
		else
		% if the last vertex is not on the same side of the line

			if(1<(i1)) 
			% delete all vertices from 1 to (i1-1)
			%borderxtmp.erase(borderxtmp.begin(), borderxtmp.begin() + (i1));
			%borderytmp.erase(borderytmp.begin(), borderytmp.begin() + (i1));
			borderxtmp(1:i1-1) = [];
			borderytmp(1:i1-1) = [];
			end

			%if((i2+1)<=(n-1)) 
			if((i2+1)<=n) 
			% delete all vertices from (i2+1) to end
			%borderxtmp.erase(borderxtmp.begin() + (i2+1), borderxtmp.end());
			%borderytmp.erase(borderytmp.begin() + (i2+1), borderytmp.end());
			borderxtmp(i2+1:end) = [];
			borderytmp(i2+1:end) = [];
			end
		end
	else 
	% if the first vertex is on the same side as that of the line

		if(ind(n)==1)
		% the last vertex is on the same side as that of the line
		
			% loop from vertices i to n-1
			for i=1:n-1
			
				% vertex i is on the same side but vertex i+1 is on the wrong side
				if((ind(i)==1) && (ind(i+1)==0)) 
				
					% store i+1 as i1
					i1 = i+1;

					% construct line2 from vertex i to i+1
					line2.x1 = borderx(i);
					line2.y1 = bordery(i);
					line2.x2 = borderx(i+1);
					line2.y2 = bordery(i+1);
					% compute intersection
					[ix,iy] = computeintersection(line, line2);
				
					% update the boundary variables
					borderxtmp(i+1) = ix;
					borderytmp(i+1) = iy;
				end
				% vertex i is on the wrong side but vertex i+1 is on the same side
				if((ind(i)==0) && (ind(i+1)==1)) 

					% store i as i2
					i2 = i;

					% construct line2 from vertex i to vertex i+1
					line2.x1 = borderx(i);
					line2.y1 = bordery(i);
					line2.x2 = borderx(i+1);
					line2.y2 = bordery(i+1);
					% compute intersection
					[ix,iy] = computeintersection(line, line2);
				
					% update the boundary variables
					if(i1==i2) 
						%borderxtmp.push_back(ix);
						%borderytmp.push_back(iy);
						%borderxtmp.insert(borderxtmp.begin()+(i+1),ix);
						%borderytmp.insert(borderytmp.begin()+(i+1),iy);
						borderxtmp = [borderxtmp(1:i)  ix  borderxtmp(i+1:end)];
						borderytmp = [borderytmp(1:i)  iy  borderytmp(i+1:end)];
					else 
						borderxtmp(i) = ix;
						borderytmp(i) = iy;
					end
				end
			end

			% remove elements from i1 to i2
			if((i1+1)<(i2)) 
				%borderxtmp.erase(borderxtmp.begin() + (i1+1), borderxtmp.begin() + (i2));
				%borderytmp.erase(borderytmp.begin() + (i1+1), borderytmp.begin() + (i2));
				borderxtmp(i1+1:i2-1) = [];
				borderytmp(i1+1:i2-1) = [];
			end

		else
		% the last vertex is on the wrong side as that of the line

			for i=1:n-1

				% vertex i is on the same side of the line but the vertex i+1 is on the wrong side of the line
				if((ind(i)==1) && (ind(i+1)==0)) 

					% store i+1 as i1
					i1 = i+1;

					line2.x1 = borderx(i);
					line2.y1 = bordery(i);
					line2.x2 = borderx(i+1);
					line2.y2 = bordery(i+1);
					% compute intersection
					[ix,iy] = computeintersection(line, line2);
				
					% update the boundary variables
					borderxtmp(i+1) = ix;
					borderytmp(i+1) = iy;
				end
			end

			i2 = n;

			% compute intersection
			line2.x1 = borderx(1);
			line2.y1 = bordery(1);
			line2.x2 = borderx(n);
			line2.y2 = bordery(n);
			[ix,iy] = computeintersection(line, line2);
				
			% update the boundary variables
			if(i1==i2) 
				%borderxtmp.push_back(ix);
				%borderytmp.push_back(iy);
				borderxtmp(end+1) = ix;
				borderytmp(end+1) = iy;
			else 
				%borderxtmp[n-1] = ix;
				%borderytmp[n-1] = iy;
				borderxtmp(n) = ix;
				borderytmp(n) = iy;

				% remove elements from i1 to i2
				if((i1+1)<(i2))
					%borderxtmp.erase(borderxtmp.begin()+i1+1,borderxtmp.begin()+i2);
					%borderytmp.erase(borderytmp.begin()+i1+1,borderytmp.begin()+i2);
					borderxtmp(i1+1:i2-1) = [];
					borderytmp(i1+1:i2-1) = [];
				end
			end
		

		end
	end
	borderx = borderxtmp;
	bordery = borderytmp;

end

%% find intersection of two lines assuming they intersect ie. they are not parallel
function [ix, iy] = computeintersection(line1, line2)

	TOL = 0.0000001;

	condition1 = 0;
	condition2 = 0;

	% check if lines have slope infinity
	if(((line1.x2 - line1.x1)<TOL) && ((line1.x2 - line1.x1)>-TOL))
		condition1 = 1;
	end
	if(((line2.x2 - line2.x1)<TOL) && ((line2.x2 - line2.x1)>-TOL))
		condition2 = 1;
	end

	if(condition1==1 || condition2==1)
	
		% both conditions will not be 1 at the same time since it is assumed that the lines are not parallel
		if(condition1==1)
			m2 = (line2.y2 - line2.y1)/(line2.x2 - line2.x1);
			ix = line1.x1;
			iy = line2.y1 + m2*(ix-line2.x1);
		else % condition2 = 1
			m1 = (line1.y2 - line1.y1)/(line1.x2 - line1.x1);
			ix = line2.x1;
			iy = line1.y1 + m1*(ix-line1.x1);
		end
	else
		% slopes of the two lines
		m1 = (line1.y2 - line1.y1)/(line1.x2 - line1.x1);
		m2 = (line2.y2 - line2.y1)/(line2.x2 - line2.x1);

		% computing the intersection
		ix = (m2*line2.x1 - m1*line1.x1 + line1.y1 - line2.y1)/(m2-m1);
		iy = line1.y1 + m1*(ix-line1.x1);
	end

end
