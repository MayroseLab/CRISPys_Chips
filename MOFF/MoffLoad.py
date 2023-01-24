"""MOFF imports and model loading"""

from keras import models
from json import loads
import globals
# import keras
# import tensorflow.python
from typing import List
from MOFF.MOFF_prediction import MOFF_score


CODE_PATH = globals.CODE_PATH
def load_moff():
    """
    Function to load MOFF algorithm
    """
    globals.moff_mtx1 = loads(open(f"{CODE_PATH}/MOFF/StaticFiles/M1_matrix_dic_D9").read())
    globals.moff_mtx2 = loads(open(f"{CODE_PATH}/MOFF/StaticFiles/M2_matrix_smooth_MLE").read())
    globals.moff_loaded_model = models.load_model(rf"{CODE_PATH}/MOFF/StaticFiles/GOP_model_3.h5")
    return moff_modified

def moff_modified(candidate_lst: List, target_lst: List) -> List[float]:
    """
    Calling MOFF algorithm

    :param candidate_lst:
    :param target_lst:
    :return:
    """
    scores = MOFF_score(globals.moff_mtx1, globals.moff_mtx2, candidate_lst, target_lst)
    return [1 - score for score in scores]
