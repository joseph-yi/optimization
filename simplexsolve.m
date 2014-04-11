function [x objective c ys] = simplexsolve(A,b,c)
%
%

% input dimensions
%

% setup LP
%flg = 0: convergence
%flg = 1: LP is unbounded below
%flg = 2: LP degenerate


flg       = 0;
[m,n]     = size(A);
mask	  = find(b<0);
A(mask,:) = -A(mask,:);
b(mask)   = -b(mask);
%need a while -- for phase 1...
%i need to run phase i, and if basical optimal is reached,
%pass onto phase ii (more modifications to raise errors etc needed)
%label everything (variable)p1/p2 for variable phase1

Ap1         =[A eye(m)];
[mp1, np1 ] = size(Ap1);
xp1         = [zeros(n,1);b];
Bp1         = find(xp1>0);
bp1         = b;
cp1         = [zeros(n,1); ones(m,1)];
Tp1         =Ap1(:,Bp1)\[bp1 eye(m)];
yp1         =Tp1(:,2:end)'*cp1(Bp1);
Tp2         = [Tp1; [cp1'*xp1,yp1']];


%in professor's code, The tp1 and Tp2 were just one variable T
%we have initialized a revised simplex tableau
%take note that we need to use Tp2 NOT Tp1

%Starting Simplex PHASE I

f = ['Starting Phase I Simplex...'];
disp(f);
disp('Initial Basis is');
disp(Bp1');
disp('Initial solution xp1, and cp1-Ap1^t*yp1 and their product');
disp([xp1      cp1-Ap1'*yp1         xp1.*(cp1-Ap1'*yp1)]);

simplex_p1   = 1;
ITERp1      = 0;
objp1       = cp1'*xp1;
Bp1OLD = Bp1;


%determine next s and r values)
while (simplex_p1 ==1)
	%in revised tableau, note that first column is for b, and rest
	%is for the y, yp2 takes the bottom values of table after
	%first index
	yp1        = Tp2(end,2:end)';
	%find the smallest value of in cp1-Ap1'*y
	%first index gives value of z min, second index gives s=column
	[zminp1, sp1] = min(cp1-Ap1'*yp1);
	%we need values of t_(1,n
	tp1       = Tp2(1:end-1,2:end)*Ap1(:,sp1)
	rp1 = Revisedgetr(np1,sp1,Bp1,Tp2,tp1)
	
    
	[Tp2, B1p1, flgp1] = RevisedSimplexTableau(Bp1,rp1,sp1,tp1,zminp1,Tp2)
	
	

	%%question here
	ITERp1=ITERp1 + 1;
	%xp1 = [zeros(n,1);b]
    xp1         = [zeros(n+m,1)];
	%%above was to initialize vector
	%% now we fill in answers for xp1 (first columnb of tableau)
	xp1(Bp1) = Tp2(1:end-1,1);
	
	obj1_p1 = cp1'*xp1;
	%checking for convergence
	%following function gives whether or not a new basis is introduced
	%by comparing which values are not in both B1p1, 
	
	%%%need to make sure none of the basis = 
	%%%DOING THE WRONG THING
	
	FF_p1 = setdiff(B1p1,Bp1);
	%if length(ff_p1)==0
	if (length(FF_p1)==0| (abs(obj1_p1)<1e-14 & ITERp1 > 1))
		disp('Simplex Phase I has converged')
		
		
		simplex_p1 = 0;
		disp('Displaying Optimal Basis');
		
		disp(Bp1');
		pause(5);
		disp('Initial solution xp1, and cp1-Ap1^t*yp1 and their product');
		disp([xp1      cp1-Ap1'*yp1         xp1.*(cp1-Ap1'*yp1)]);
	    continue;
	end
	Bp1 = B1p1;
	objp1 = obj1_p1;
	disp('Current Basis is');
	disp(Bp1');
	%pause(1);
end
simplex=1
if cp1'*xp1>0,

	disp('LP is degenerate - no feasible solution exists');
	disp(xp1)
	disp(cp1)
	cp1'*xp1
	simplex=0
	end



%degenerate check

%DC = ismember(Bp1OLD,Bp1);

%if sum(DC) >0,
		
	%	disp('LP is degenerate: no feasible solution exists');
		%report optimal solution and obj value
		
	%	disp('Displaying optimal solution x, obj value, and B');
	%	disp([xp1]);
	%	disp([cp1'*xp1]);
	%	disp(Bp1);
	%	
 %       simplex = 0;
		
		
%end


[m n] = size(A);
x=xp1(1:n);


B=find(x>0);
T = A(:,B)\[b eye(m)];
y = T(:,2:end)'*c(B);
size(y)
size(T)
size(c'*x)
T = [T;[c'*x,y']];
%
% Starting Simplex Method Phase II



if simplex~=0,

f = ['Starting Phase II Simplex Iteration... '];
format short g;
disp(f);
disp('Initial Basis is');
disp(B');
disp('Displaying Initial solution x, c-A^T*y and their componentwise product');
disp([x c-A'*y x.*(c-A'*y)]);
simplex = 1;
ITER = 0;
obj = c'*x;
%pause(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%QUESTION

%%%%%%%%%%%%%%%%%%%%%%%%%QUESTION

while (simplex == 1)
%
% determine the next s and r values.
%
   y        = T(end,2:end)';
   [zmin,s] = min(c-A'*y); 
   t        = T(1:end-1,2:end)*A(:,s);
   r = Revisedgetr(n,s,B,T,t);
   
   if (r < 1)
       disp('LP has no lower bound');
	   %A(:,B)\A(:,s);
	   fprintf('Take solution and add\n')
	   disp(A(:,s))
	   
	   %disp('Take solution and add column %d', A(:,s))
       simplex = 0;
       continue;
   end
   x   = zeros(n,1);
   x(B)= T(1:end-1,1);
   ITER = ITER + 1;
   f = ['Iteration ', num2str(ITER), ' Obj ', num2str(c'*x), '. Smallest component in c-A^T*y: ', ... 
         num2str(zmin), ' at s =', num2str(s), '. Component r = ', num2str(r), ' out of basis'];
   disp(f);
   obj1 = c'*x;
%
% update the revised simplex tableau.
%
   [T,B1,flg]=RevisedSimplexTableau(B,r,s,t,zmin,T); 
%%do i need this? degeneracy should be taken care of in phase 1 check
   if (flg == 1)
       disp('LP is degenerate');
       simplex = 0;
       continue;
   end
%
% check for convergence.
%
   FF  = setdiff(B1,B);
   if (length(FF) == 0 | (abs(obj1-obj) < 1e-14 & ITER  > 1))
       disp('Simplex Method has converged');
	   disp(abs(obj1-obj));
	   disp(setdiff(B1,B));
	   disp('hello i am a cat');
       simplex = 0;
       disp('Displaying Optimal Basis');
       disp(B');
       disp('Displaying Optimal solution x, c-A^T*y and their componentwise product');
       disp([x c-A'*y x.*(c-A'*y)]);
	   
       continue;
   end
   B   = B1;
   obj = obj1;
   disp('Current Basis is');
   disp(B')
   %pause(1);
end
objective=c'*x;
ys =T(end,2:end)';


