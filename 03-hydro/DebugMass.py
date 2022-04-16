import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio
from welib.yams.windturbine import FASTWindTurbine



fstFilename = 'TS_MD0HD1SD0_F101010/Main.fst'; 

WT = FASTWindTurbine(fstFilename, twrShapes=[0,2], algo='OpenFAST')
# print(WT.twr)
print('>>>>',np.around(WT.twr.mass     ,3), 256128.857)
print('>>>>',np.around(WT.bld[0].mass  ,3), 17730.328)
print('>>>>',np.around(WT.rot.mass     ,3), 96193.983)
print('>>>>',np.around(WT.fnd.mass     ,3), 5123033.000)
print('>>>>',np.around(WT.WT_rigid.mass,3), 5584057.839)
print('>>>>',np.around(WT.WT_rigid.mass,3), 5584059.8)

if __name__ == '__main__':
    pass
