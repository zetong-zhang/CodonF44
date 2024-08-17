# Author: Zhang Zetong
# Date: 2024-08-17
import pickle

with open("CodonF44_Model.pkl", "rb") as file:
    model = pickle.load(file)


def SVM(vec):
    result = model.predict([vec])
    return result[0]
