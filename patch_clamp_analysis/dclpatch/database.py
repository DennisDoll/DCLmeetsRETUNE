from pathlib import Path
from abc import ABC, abstractmethod
from pandas import DataFrame

import pandas as pd
import numpy as np


class Database:
    
    """
    The database is supposed to hold all information about all recorded cells that were added to it.
    These information exclude the raw data (only contain filepaths to raw data),
    but include general metadata like celltype, stimulus type, brain region, pharamacological treatment.
    Potentially, some intermediate data could be added (like how the cell reacted upon stimulation, see dashboard).
    Overall, the database should allow the selection of cells based on mix-and-match criteria 
    for further (statistical) analyses and plotting.
    """
    
    def __init__(self, root_dir: Path, raw_data_dir: Path):
        self.root_dir = root_dir
        self.subdirectories = Subdirectories(root_dir = root_dir)
        self.global_cell_id = 0 

    def load_database_from_disk(self):
        # Option to load previously created and saved database
        pass
    
    def save_database_to_disk(self):
        # Save all information to disk
        pass
    
    def add_new_cell_recording(self, cell_recordings_dir: Path):
        new_recording = CellRecording()
        recording_overview = new_recording.create_recordings_overview(cell_recordings_dir = cell_recordings_dir,
                                                                      global_cell_id = self.global_cell_id)
        if hasattr(self, 'recorded_cells') == False:
            self.recorded_cells = recording_overview
        else:
            self.recorded_cells = pd.concat([self.recorded_cells, overview])
        
        self.global_cell_id += 1
        
        # Trigger update of mix-and-match categories?
        
    

class CellRecording:
    
    def create_general_metadata_df(path_to_recordings_dir: Path) -> DataFrame:
        metadata_df = pd.read_excel(path_to_recordings_dir.joinpath(f'{path_to_recordings_dir.name}.xlsx'),
                                    sheet_name = 'General information')
        general_metadata = {'date': path_to_recordings_dir.name[:10],
                            'session_cell_id': path_to_recordings_dir.name[path_to_recordings_dir.name.rfind('_') + 1:],
                            'mouse_line': metadata_df['Animal line'][0],
                            'brain_region': metadata_df['Region'][0],
                            'cell_type': metadata_df['Type'][0],
                            'sex': metadata_df['Sex'][0]}
        return pd.DataFrame(general_metadata, index=[0])
    
    
    def create_stimulation_paradigms_df(path_to_recordings_dir: Path) -> DataFrame:
        filepaths_stimulation_recordings = list()
        for elem in path_to_recordings_dir.iterdir():
            if 'datapoints' not in elem.name:
                if elem.name.endswith('.csv'):
                    filepaths_stimulation_recordings.append(elem)
        stimulation_paradigms = {'stimulation_string': list(),
                                 'stimulation_frequency-Hz': list(),
                                 'stimulation_duration-ms': list(), 
                                 'filepath_detected_events': list()}
        for filepath in filepaths_stimulation_recordings:
            filename = filepath.name
            filename = filename.replace('.csv', '')
            yyyy_mm_dd_ = filename[:11]
            filename = filename.replace(yyyy_mm_dd_, '')
            cell_id = filename[:filename.index('_')]
            filename = filename.replace(cell_id + '_', '')
            if 'Hz' in filename:
                stimulation_frequency = int(filename[:filename.index('Hz')])
                stimulation_duration = filename[filename.index('Hz') + 2 :]
                if 'ms' not in stimulation_duration:
                    stimulation_duration = int(stimulation_duration[:stimulation_duration.find('s')]) * 1000
                else:
                    stimulation_duration = int(stimulation_duration[:stimulation_duration.find('ms')])
                stimulation_paradigm = f'{stimulation_frequency}-Hz_for_{stimulation_duration}-ms'
            elif 'Bsl' in filename:
                stimulation_frequency = np.NaN
                stimulation_duration = np.NaN
                stimulation_paradigm = 'baseline'
            else:
                print(f'Warning: stimulation paradigm could not be identified for: {filepath.name}')
                stimulation_frequency = np.NaN
                stimulation_duration = np.NaN
                stimulation_paradigm = 'unknown'
            stimulation_paradigms['stimulation_string'].append(stimulation_paradigm)
            stimulation_paradigms['stimulation_frequency-Hz'].append(stimulation_frequency)
            stimulation_paradigms['stimulation_duration-ms'].append(stimulation_duration)
            stimulation_paradigms['filepath_detected_events'].append(filepath)
        return pd.DataFrame(stimulation_paradigms)
    
    
    def create_recordings_overview(self, cell_recordings_dir: Path, global_cell_id: int) -> DataFrame:
        general_metadata = self.create_general_metadata_df(cell_recordings_dir)
        stimulation_paradigms = self.create_stimulation_paradigms_df(cell_recordings_dir)
        stimulations_adjusted_general_metadata = pd.concat([general_metadata]*stimulation_paradigms.shape[0], ignore_index=True)
        recordings_overview = pd.concat([stimulations_adjusted_general_metadata, stimulation_paradigms], axis=1)
        recordings_overview['global_cell_id'] = str(global_feature_id).zfill(4)
        return recordings_overview


class Subdirectories:
    
    def __init__(self, root_dir: Path):
        self.root_dir = root_dir
        self.create_missing_subdirectories()
        self.assign_subdirectories_as_attributes()
        
    def create_missing_subdirectories(self):
        # check for each element in a list of subdirs, whether they exist --> create if not
        pass
    
    def assign_subdirectories_as_attributes(self):
        # use list of subdirs and set the path to each as attribute
        pass