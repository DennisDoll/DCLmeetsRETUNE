from pathlib import Path

from .database import Database
from .analysis import SingleCellAnalysis
from typing import List, Dict



class PatchProject:
    
    def __init__(self, root_dir: Path):
        self.database = Database(root_dir = root_dir)
    
    def add_cell_to_database(self, path_to_cell_recordings_dir: Path, overwrite: bool=False):
        self.database.add_new_cell_recording(cell_recordings_dir = path_to_cell_recordings_dir, overwrite = overwrite)
        
    def compare_on_single_cell_level(self, global_cell_id: str, show: bool=True, save: bool=False):
        single_cell_analysis = SingleCellAnalysis(database = self.database)
        single_cell_analysis.run_analysis(global_cell_id = global_cell_id, show = show, save = save)
    
    def get_global_cell_ids_matching_criteria(self, criteria: Dict) -> List:
        pass
    
    def compare_on_group_level(self, global_cell_ids_group_a: List, global_cell_ids_group_b: List):
        pass