import os
from .database import Database
from typing import List, Dict, Tuple, Optional

"""

Longterm: should the other "main" classes like "preprocessor" or "segmentor" be integrated here in main.py?

"""

class Project:
    def __init__(self, user_input: Dict):
        self.project_root_dir = user_input['project_root_dir']
        self.database = Database(user_input)
        
    def save_status(self) -> None:
        self.database.save_all()
    
    def load_status(self) -> None:
        self.database.load_all()
    
    def preprocess(self, file_ids: Optional[List]=None, overwrite: bool=False) -> None:
        from .preprocessing import Preprocessor
        
        if file_ids == None:
            file_ids = self.database.get_file_ids_to_process(process_tracker_key = 'preprocessing_completed', overwrite = overwrite)
        
        for file_id in file_ids:
            preprocessing_object = PreprocessingObject(database = self.database, file_id = file_id)
            preprocessing_object.run_all_preprocessing_steps()
            preprocessing_object.save_preprocessed_images_on_disk()
            preprocessing_object.save_preprocessed_rois_in_database()
            #update database about the progress
            del preprocessing_object

                
        preprocessor = Preprocessor(file_ids, self.database)
        self.database = preprocessor.run_individually()
    
    def run_segmentation(self, file_ids: Optional[List]=None) -> None:
        from .segmentation import Segmentor

        segmentor = Segmentor(self.database, file_ids)
        self.database = segmentor.run_all()

    def run_quantifications(self, file_ids: Optional[List]=None) -> None:
        from .quantifications import Quantifier
        
        if 'quantification_completed' not in self.database.file_infos.keys():
            self.database.add_new_key_to_file_infos('quantification_completed')
        if file_ids == None:
            all_file_ids = self.database.file_infos['file_id']
            quantification_status = self.database.file_infos['quantification_completed']
            file_ids = [elem[0] for elem in zip(all_file_ids, quantification_status) if elem[1] == False or elem[1] == None]
        quantifier = Quantifier(self.database, file_ids)
        self.database = quantifier.run_all()
        
    def run_inspection(self, file_id: str, inspection_strategy):
        from .inspection import InspectionStrategy
        
        inspection_strategy.run(self.database, file_id)
        
    def remove_file_id_from_project(self, file_id: str):
        self.database.remove_file_id_from_project(file_id = file_id)

        
        
      
        

        
        
        
        
        
        
        


            
            
            
