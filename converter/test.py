import json

with open("config.json") as f:
    cfg = json.load(f)

import subprocess
for cmd in cfg["commands"]:
    subprocess.run(cmd, shell=True)
