from pathlib import Path

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from .database import Database


class SingleCellAnalysis:
    
    def __init__(self, database: Database):
        self.database = database
        
    def run_analysis(self, global_cell_id: str, show: bool, save: bool):
        