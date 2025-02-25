import pytest
import numpy as np
import math

import sys
from pathlib import Path
root = Path(__file__).parent.parent
sys.path.append(str(root))

from scripts.experiments.BFQ_Simulation import main


@pytest.mark.benchmark
def test_main_function(benchmark):
    benchmark(main)
