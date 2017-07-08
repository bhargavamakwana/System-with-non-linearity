%please read the code before executing it. Refer to the following paper in case if
%needed "global stability of relay feedback systems" by Jorge M.goncalves,
%IEEE transactions on automatic control, VOL 46, NO 4 April 2001.

clc;
clear all;
t=0;
k=0:0.001:10;
% hystresis parameter 
d=0.25;
%numerator for transfer function
b=[1 -5 -11 -5 250 250];
%denominator for transfer function
a=[1 12 57 140 195 150 50];
%from transfer function to state space
[A,B,C,D]=tf2ss(b,a);
%identity matrix
I=eye(6);
t_star=2.094;

% comment x* to fnd t* from the graph of g. And set t* for further calculations.
%uncomment for other observations.
x_s=inv(expm(A*t_star)+I)*(expm(A*t_star)-I)*inv(A)*B;
for i=1:10001
    ea=expm(t*A);  % matrix exponential 
    
    % limit cycle funtion zero at t*
    % find first the t* time to know the period. As it is used in further
    % calculations
    %for reference for y and g refer equation 4 of the paper
    
    g=C*(inv(ea+I))*(ea-I)*(inv(A))*B-d;   % g(t)
    G(:,i)=g;
    %comment from here 
    y=C*(ea*((x_s)-(inv(A)*B))+(inv(A)*B))+d; % output as y(t)+d 
    Y(:,i)=y;
    
    %for vt refer theorem 4.1 of the paper.
    vt=(expm(t*A)*x_s)-ea*(inv(A)*B)-(expm(A*t_star)*x_s)+(expm(A*t_star)*(inv(A)*B));
    
    %Jacobian of the poincare maps computed at every iteration. 
    H=((inv(C*vt)*(vt*C))-I)*ea;
    
    %solve the LMI equations
    
    setlmis([])
    p = lmivar(1,[6 1]);
    lmiterm([-1 1 1 p],1,1);
    lmiterm([1 1 1 p],H',H);
    lmiterm([2 1 1 p],1,1);
    lmis = getlmis;
    [tmin,xfeas] = feasp(lmis);
    
    %get the P matrix.
    P = dec2mat(lmis,xfeas,p);
    deriv=P-(H'*P*H);
    
    %find the eigen values to validate the inequality 
    eig_values=eig(deriv);
    
    % if min is positive then stable.
    value=min(double(eig_values));
    V(:,i)=value;
    %comment till here
    
    t=t+0.001;
end
subplot(3,1,1);

%plot g and find t* for g=0.
plot(k,G,'r');
subplot(3,1,2);

%plot y.
plot(k,Y,'y');
subplot(3,1,3);

%plot of the minimum eigen values 
plot(k,V,'g');