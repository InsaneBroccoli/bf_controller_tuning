function out = fix_signal(sig)

out = sig;
idx = find(movsum(sig == 0, 500) == 500, 1);
out(idx:end) = 0;
idx = find(out == 0, 1);
out(1:idx) = 0;

end

