import sys
import os
import pytest

# Add the current directory to sys.path
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

if __name__ == "__main__":
    # Run the tests
    pytest.main(["tests", "-v"])
