% run_json_commands.m  (script or put inside a function if you prefer)
cfg = jsondecode(fileread("config.json"));

cmds = cfg.commands;
if isstring(cmds), cmds = cellstr(cmds); end
if ischar(cmds),   cmds = {cmds};        end

for k = 1:numel(cmds)
    % system(...) uses the OS shell:
    % - Windows: cmd.exe
    % - macOS/Linux: /bin/sh
    system(cmds{k});   % prints output to Command Window; no error thrown on nonzero exit
end