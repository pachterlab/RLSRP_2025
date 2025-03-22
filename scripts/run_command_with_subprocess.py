import subprocess
import argparse
import ast

parser = argparse.ArgumentParser()
parser.add_argument("items", type=str, help="A list passed as a string")
args = parser.parse_args()

try:
    command = ast.literal_eval(args.items)
    subprocess.run(command, check=True)
except Exception as e:
    command = args.items
    raise e  # subprocess.run(command, check=True, shell=True)  # keep commented out so that I don't support string commands
