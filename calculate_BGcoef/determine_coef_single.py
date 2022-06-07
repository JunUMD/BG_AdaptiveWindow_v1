# -*- coding: utf-8 -*-
"""                 
Created on 12/25/2020
                    
@author: Jun@UMD             
"""


def determine_coef_single(input1):
    import calculate_BGcoef.determine_coef_pythonlib.determine_coef_class

    [nedt,ifr,numxgrd,numygrd,
     win_file,uvs_file,
     noise_lev,save_file]=input1
 

    determine_coef_obj=calculate_BGcoef.determine_coef_pythonlib.determine_coef_class. \
                       determine_coef_class(nedt,ifr,numxgrd,numygrd,
                                            win_file,uvs_file,
                                            noise_lev,save_file)
    determine_coef_obj.ingest()
    
    determine_coef_obj.process()

