function output = expma(x)
  A = logical(x>300);
  B = logical(x<(-300));
  output = exp(x);
  output(A) = exp(300);
  output(B) = exp(-300);
end

    