import importlib_resources
from io import StringIO
import pandas as pd


def import_package_file(pkg_name, filename):
    ref = importlib_resources.files(pkg_name) / filename
    with importlib_resources.as_file(ref) as path:
        with open(path, 'r') as f:
            content = f.read()

    # convert a csv (as text) into a pandas csv
    if filename[-4:] == ".csv":
        content = pd.read_csv(StringIO(content), sep=",")
    
    return content

    