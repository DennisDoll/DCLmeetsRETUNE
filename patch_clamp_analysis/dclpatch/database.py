from pathlib import Path
from abc import ABC, abstractmethod
from typing import List
from pandas import DataFrame
import pickle
import os
from datetime import datetime

import pandas as pd
import numpy as np

SORTED_COLUMNS = ['global_cell_id', 'date', 'session_cell_id', 'mouse_line', 'sex', 'brain_region', 'cell_type', 
                  'stimulation_string', 'stimulation_frequency-Hz', 'stimulation_duration-ms', 'filepath_detected_events']


def listdir_nohidden(path: Path) -> List:
    return [f for f in os.listdir(path.as_posix()) if f.startswith('.') == False]

class Database:
    
    """
    The database is supposed to hold all information about all recorded cells that were added to it.
    These information exclude the raw data (only contain filepaths to raw data),
    but include general metadata like celltype, stimulus type, brain region, pharamacological treatment.
    Potentially, some intermediate data could be added (like how the cell reacted upon stimulation, see dashboard).
    Overall, the database should allow the selection of cells based on mix-and-match criteria 
    for further (statistical) analyses and plotting.
    """
    
    def __init__(self, root_dir: Path):
        self.root_dir = root_dir
        self.subdirectories = Subdirectories(root_dir = root_dir)
        self.global_cell_id = 0
        self.cells_in_database = list()

    def load_database_from_disk(self):
        # Option to load previously created and saved database
        pass
    
    def save_database_to_disk(self):
        # Save all information to disk
        pass


    def list_all_column_values(self, global_cell_id: str, column_name: str) -> List:
        if hasattr(self, 'cell_recordings_metadata'):
            if global_cell_id in self.cell_recordings_metadata['global_cell_id'].unique():
                values = list(self.cell_recordings_metadata.loc[self.cell_recordings_metadata['global_cell_id'] == global_cell_id, column_name].values)
                return values
            else:
                raise ValueError(f'{global_cell_id} is not in the database yet!')
        else:
            print('No entries in the database yet')
                                 
    
    def add_new_cell_recording(self, cell_recordings_dir: Path, overwrite: bool):
        increase_global_cell_id_count = True
        if overwrite:
            if cell_recordings_dir in self.cells_in_database:
                global_cell_id = self.cell_recordings_metadata.loc[self.cell_recordings_metadata['filepath_detected_events'].str.startswith(cell_recordings_dir.as_posix()), 'global_cell_id'].iloc[0]
                increase_global_cell_id_count = False
                self.cell_recordings_metadata.drop(self.cell_recordings_metadata.loc[self.cell_recordings_metadata['global_cell_id'] == global_cell_id].index, inplace = True)
                self.cells_in_database.remove(cell_recordings_dir)
            else:
                global_cell_id = self.global_cell_id
                

        if cell_recordings_dir not in self.cells_in_database:
            new_recording = CellRecording()
            if overwrite == False:
                global_cell_id = self.global_cell_id
            recording_overview = new_recording.create_recordings_overview(cell_recordings_dir = cell_recordings_dir,
                                                                          global_cell_id = global_cell_id)
            if hasattr(self, 'cell_recordings_metadata') == False:
                for column_name in list(recording_overview.columns):
                    if column_name not in SORTED_COLUMNS:
                        print(f'Warning: {column_name} not defined in SORTED_COLUMNS.')
                self.cell_recordings_metadata = recording_overview[SORTED_COLUMNS]
            else:
                for column_name in list(recording_overview.columns):
                    if column_name not in list(self.cell_recordings_metadata.columns):
                        print(f'Warning: {column_name} represents a new column. Consider updating of all previously added cell recordings.')
                self.cell_recordings_metadata = pd.concat([self.cell_recordings_metadata, recording_overview], ignore_index = True)
            self.cell_recordings_metadata.sort_values(by=['global_cell_id'], inplace = True, ignore_index = True)
            self.cells_in_database.append(cell_recordings_dir)
            if increase_global_cell_id_count:
                self.global_cell_id += 1
        
            # Trigger update of mix-and-match categories?
        else:
            note_line1 = f'Note: The recordings in "{cell_recordings_dir.as_posix()}" were already added to the database.\n'
            note_line2 =  '      Consider passing "overwrite = True" when calling this function in order to update the information.'
            print(note_line1 + note_line2)
            
            
    def save_all(self):
        self.save_all_attributes_as_pickle()
        self.save_cell_recordings_metadata_as_csv()

    def save_all_attributes_as_pickle(self):
        project_summary = self.__dict__.copy()
        filepath = f'{self.subdirectories.project_summary.as_posix()}/{datetime.now().strftime("%Y_%m_%d")}_dclpatch_project_summary.p'        
        with open(filepath, 'wb') as io:
            pickle.dump(project_summary, io)

    def save_cell_recordings_metadata_as_csv(self):
        filepath = f'{self.subdirectories.project_summary.as_posix()}/{datetime.now().strftime("%Y_%m_%d")}_dclpatch_database_overview.csv'  
        self.cell_recordings_metadata.to_csv(filepath)
    
    
    def load_all(self):
        result_files = [fname for fname in listdir_nohidden(self.subdirectories.project_summary) if fname.endswith('dclpatch_project_summary.p')]
        result_files.sort(reverse = True)
        with open(f'{self.subdirectories.project_summary.as_posix()}/{result_files[0]}', 'rb') as io:
            project_summary = pickle.load(io)

        for key, value in project_summary.items():
            if hasattr(self, key) == False:
                setattr(self, key, value)
        
    

class CellRecording:
    
    def create_general_metadata_df(self, path_to_recordings_dir: Path) -> DataFrame:
        metadata_df = pd.read_excel(path_to_recordings_dir.joinpath(f'{path_to_recordings_dir.name}.xlsx'),
                                    sheet_name = 'General information')
        general_metadata = {'date': path_to_recordings_dir.name[:10],
                            'session_cell_id': path_to_recordings_dir.name[path_to_recordings_dir.name.rfind('_') + 1:],
                            'mouse_line': metadata_df['Animal line'][0],
                            'brain_region': metadata_df['Region'][0],
                            'cell_type': metadata_df['Type'][0],
                            'sex': metadata_df['Sex'][0]}
        return pd.DataFrame(general_metadata, index=[0])
    
    
    def create_stimulation_paradigms_df(self, path_to_recordings_dir: Path) -> DataFrame:
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
            stimulation_paradigms['filepath_detected_events'].append(filepath.as_posix())
        return pd.DataFrame(stimulation_paradigms)
    
    
    def create_recordings_overview(self, cell_recordings_dir: Path, global_cell_id: int) -> DataFrame:
        general_metadata = self.create_general_metadata_df(cell_recordings_dir)
        stimulation_paradigms = self.create_stimulation_paradigms_df(cell_recordings_dir)
        stimulations_adjusted_general_metadata = pd.concat([general_metadata]*stimulation_paradigms.shape[0], ignore_index=True)
        recordings_overview = pd.concat([stimulations_adjusted_general_metadata, stimulation_paradigms], axis=1)
        recordings_overview['global_cell_id'] = str(global_cell_id).zfill(4)
        return recordings_overview


class Subdirectories:
    
    def __init__(self, root_dir: Path):
        self.root_dir = root_dir
        self.create_missing_subdirectories()
        
    def create_missing_subdirectories(self):
        subdirectories = listdir_nohidden(self.root_dir)
        
        try: self.single_cell_analyses = self.root_dir.joinpath([elem for elem in subdirectories if 'single_cell' in elem][0])
        except:
            self.single_cell_analyses = self.root_dir.joinpath('single_cell_analysis')
            os.mkdir(self.single_cell_analyses.as_posix()) 
        
        try: self.group_analyses = self.root_dir.joinpath([elem for elem in subdirectories if 'group' in elem][0])
        except:
            self.group_analyses = self.root_dir.joinpath('group_analysis')
            os.mkdir(self.group_analyses.as_posix()) 
            
        try: self.project_summary = self.root_dir.joinpath([elem for elem in subdirectories if 'project_summary' in elem][0])
        except:
            self.project_summary = self.root_dir.joinpath('project_summary')
            os.mkdir(self.project_summary.as_posix()) 