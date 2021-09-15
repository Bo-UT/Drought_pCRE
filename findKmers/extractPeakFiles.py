import os
from glob import glob
# %%
for i in os.walk('peaks'):
    print('first')
    print(i)

# %%
files = [os.path.join(dp, f) for dp, dn, filenames in os.walk('peaks/')\
     for f in filenames if os.path.splitext(f)[1] == '.narrowPeak']
dapFiles = [file for file in files if '_col_' in file.split('/')[2]]
# %%
# move and rename files to a new dierctory
for file in dapFiles:
    os.rename(file, 'DAPnarrowPeaks/{}.narrowPeak'.\
        format(file.split('/')[2].split('_')[0]))
# %%

# %%
