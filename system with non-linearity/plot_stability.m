clc;
clear all;
% write numerator and denominator coefficients for the function.
num=input('enter for numerator: ');
den=input('eneter for denominator: ');
%convert to transfer function
H=tf(num,den);
%from transfer function to satate space
[A,B,C,D]=tf2ss(num,den);
t=0:0.01:10;


d=input('enter the size of relay i.e in form of d where 2d is width: ');
nyquist(H);
plot(t,Y,'r');
hold on;
plot(t,U,'g');

