import h5py
import numpy as np
from pathlib import Path
import os
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Type, Union, cast
import yaml
from yaml import Loader, Dumper, SafeLoader
import re
import yaml
import sys
import ast

lib_to_import=sys.argv[0].replace("/scripts/create_yaml.py", "/")
sys.path.append(lib_to_import+"common")
from utils import *

config=sys.argv[1]
config=ast.literal_eval(config)

print("make yaml")
make_yaml(config)