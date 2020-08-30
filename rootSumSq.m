function [Y] = rootSumSq(X,d)
  if nargin<2
    d = 1;
  end
  Y = sqrt(sum(abs(X).^2,d));
end