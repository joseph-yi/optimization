
file_id=fopen('infile.txt');

sizes = str2num(fgetl(file_id));

n=sizes(1);
m=sizes(2);

A_matrix =	zeros(n,m);
b        =  zeros(n,1);
for i=1:(n*m),
	line        	= fgetl(file_id);
	line_vector 	= str2num(line);
	if size(line_vector)~=3,
		err=MException('there should be 3 values!')
		throw(err)
		end
	
	i               = line_vector(1);
	j               = line_vector(2);
	A(i,j)          = line_vector(3);
	end
	
for i=1:n,
	line 		= fgetl(file_id);
	line_vector = str2num(line);
	
	if size(line_vector) ~=2,
		err=MException('there should be 2 values!')
		throw(err)
		end
	b(line_vector(1)) = line_vector(2);
	
	end

last = fgetl(file_id);
if last ~= -1;
	err=MException('why is there a value here.')
	throw(err)
	end


%choose initial basis as first m components
t	=zeros(n,1);
B	=1:m;
Bc	=m+1:n;

%give initial solution
x		= A(B,:)^-1*b(B);
%give initial y
h		= A*x-b;
y(Bc)	=sign(h(Bc));
y(B)	= -(A(B,:)^(-1))'*A(Bc,:)'*y(Bc)';
%create initial tableau

if rank(A(B,:))~= min(size(A(B,:)))
	disp('matrix not invertible: degenerate case');
	break,
	end;
	
M		= A(B,:);
Minv	= A(B,:)^-1;
cost	= sum(abs(A*x-b));
Tp1 	= [x Minv; cost y(B)];
y_abs 	= abs(y(B));
%now check if we need to go to phase 2
counter = 0;
while any(abs(y(B))>1)
	counter=counter+1;
	h		=A*x-b;
	y_abs	=abs(y(B));	
	
	%get s and r
	
	all_s	= find(y_abs>1);
	s		= B(all_s(1));
	j		= find(B==s);
	e		= zeros(m,1);
	e(j)	= 1;
	%DEGENERACY TESTS
	
	if rank(A(B,:))~= min(size(A(B,:))),
	disp('matrix not invertible: degenerate case');
	break,
	end;
	
	
	%get r	
	t(Bc)		=  -sign(y(s))*(y(Bc)').*(A(Bc,:)*A(B,:)^(-1)*e);
	t(B)		= x;
	
	masking		= find(t(Bc)>0);
	[lambda r]	= min(abs(h(masking))./t(Bc(masking)));
	r 			= Bc(masking(r));
	
	%now update the basis vector
	B			= setdiff(B,s);
	B			= [B r];
	Bc 			= setdiff(1:n,B);
	%update the table
	%the u matrix
	
	u		= A(r,:)*Minv;
	%actually the transpose
	for i=1:m,
	
		if i==j,
			Minv(:,i) = Minv(:,i)/u(j); 
		else
			Minv(:,i) = Minv(:,i)-u(i)/u(j)*Minv(:,j);
		end
	end
	
	A(B,:) 	= Minv^-1
	
	x		= Minv*b(B);
	cost	= sum(abs(A*x-b));
	
	y(Bc)	= sign(h(Bc));

	y(B)	= -(Minv)'*A(Bc,:)'*y(Bc)';
	
	disp(y(B));
	disp('jam jam jam');
	disp(cost);
	%need to update y(B) at the end....
	disp('yay')
	%if doesnt converge after 10,000 steps... something probably wrong
	if (counter>10000),
		disp('i dont think my solution is converging...unknown error!')
		break
		end
	end

counter = 0
if all(abs(y(B) <= 1))
	disp('we have an optimal solution');
	disp('showing solution');
	disp(x);
	disp('showing the final basis');
	disp(B);
	disp('showing the y values');
	disp(y(B));
	disp('showing the norm');
	disp(cost);
end
file_Name 	= 'output.txt';
fid 	= fopen(file_Name,'w');
string1 = 'This is the optimal primal solution';
string2 = 'This is the optimal dual solution';
string3 = 'This is the optimal objective value';

fprintf(fid, '%s\r\n',string1);
fprintf(fid, '%f\r\n\n',x);

fprintf(fid,'%s\r\n',string2);
fprintf(fid,'%f\r\n', y(B))

fprintf(fid,'%s\r\n',string3);
fprintf(fid,'%f\r\n',cost);
fclose(fid)





