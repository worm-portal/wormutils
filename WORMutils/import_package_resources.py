import importlib_resources
from io import StringIO
import pandas as pd


def import_package_file(pkg_name, filename, import_csv_as_df=True, as_file=False):
    ref = importlib_resources.files(pkg_name) / filename

    if as_file:
        # Return the context manager for use with 'with' statement
        return importlib_resources.as_file(ref)
    else:
        # Read content as string
        with importlib_resources.as_file(ref) as path:
            with open(path, 'r') as f:
                content = f.read()
                
        # Convert CSV to pandas DataFrame if requested
        if filename[-4:] == ".csv" and import_csv_as_df:
            content = pd.read_csv(StringIO(content), sep=",")
    
    return content