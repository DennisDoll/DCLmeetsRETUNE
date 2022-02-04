from pathlib import Path
import os

from .database import Database, listdir_nohidden
from .analysis import CDFAnalysis
from typing import List, Dict



class PatchProject:
    
    def __init__(self, root_dir: Path):
        self.database = Database(root_dir = root_dir)
    
    def add_cell_to_database(self, path_to_cell_recordings_dir: Path, overwrite: bool=False):
        self.database.add_new_cell_recording(cell_recordings_dir = path_to_cell_recordings_dir, overwrite = overwrite)
        
    def remove_cell_from_database(self):
        pass
    
    def add_all_cells_to_database(self, path_with_cell_recordings_as_subdirs: Path, overwrite: bool=False):
        subdirectories = list()
        elements = listdir_nohidden(path = path_with_cell_recordings_as_subdirs)
        for elem in elements:
            subdir = f'{path_with_cell_recordings_as_subdirs.as_posix()}/{elem}'
            if os.path.isdir(subdir):
                subdirectories.append(Path(subdir))
        for cell_recordings_dir in subdirectories:
            path_to_cell_recordings_dir = Path(cell_recordings_dir)
            self.add_cell_to_database(path_to_cell_recordings_dir = path_to_cell_recordings_dir, overwrite = overwrite)
        
    def compare_on_single_cell_level(self, global_cell_id: str, show: bool=True, save: bool=False):
        single_cell_analysis = CDFAnalysis(database = self.database)
        single_cell_analysis.run_analysis(group_column = 'global_cell_id', group_id = global_cell_id, show = show, save = save)
    
    def compare_within_group(self, group_column: str, group_id: str, show: bool=True, save: bool=False):
        group_analysis = CDFAnalysis(database = self.database)
        group_analysis.run_analysis(group_column = group_column, group_id = group_id, show = show, save = save)
        
    def compare_between_groups(self):
        pass