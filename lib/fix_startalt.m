function out = fix_startalt(sig, idx)
out = sig;
out(1:idx - 1) = 0;
out(idx:end) = out(idx:end) - out(idx);
end
