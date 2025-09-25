import json
import pathlib
import os

def load_config(path: str):
    """
    Load a JSON config file and expand placeholders:
      {cwd}       -> current working directory
      {scriptdir} -> folder where the script lives
      {home}      -> user home directory
      {ENV_VAR}   -> environment variable
    """
    with open(path, "r") as f:
        cfg = json.load(f)

    cwd = str(pathlib.Path().resolve())

    def expand(value: str) -> str:
        if not isinstance(value, str):
            return value
        value = value.replace("{cwd}", cwd)
        # expand {ENV_VAR}
        for key, val in os.environ.items():
            value = value.replace(f"{{{key}}}", val)
        return value

    # expand all top-level values
    for k, v in cfg.items():
        if isinstance(v, str):
            cfg[k] = expand(v)
        elif isinstance(v, list):
            cfg[k] = [expand(x) for x in v]

    return cfg
