function out = fix_signal(sig)

out = sig;
idx = find(out == 0, 1);
out(1:idx) = 0;

end

