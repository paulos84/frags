import pandas as pd

from compounds.models import Compound

def foo():
    df = pd.DataFrame(Compound.chemical_data())
    df.to_csv('cpd_data.csv')
    return df