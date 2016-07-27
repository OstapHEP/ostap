import sys

# Check the Python version
if sys.version_info < (2, 6):
    sys.exit("Ostap requires Python 2.6 and above.")
