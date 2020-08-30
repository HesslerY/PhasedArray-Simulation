function [Y] = clamp(X,a,b)
  Y = X;
  Y(X<a)=a;
  Y(X>b)=b;
end