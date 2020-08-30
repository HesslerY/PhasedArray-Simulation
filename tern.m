function [out] = tern(sw,tr,fs)
  if sw
    out = tr;
  else
    out = fs;
  end
  if isa(out,'function_handle')
    out = out();
  end
end