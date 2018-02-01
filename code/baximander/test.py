"""
This file converts all the .mat files in current folder into readable format
"""

import scipy.io
from pprint import pprint 
import os

for subdir, dirs, files in os.walk('.'):
	for file in files:
		filepath = subdir + os.sep + file
		if not file.endswith('mat'): continue
		data = scipy.io.loadmat(file)
		filename = file.split('.')[0]
		with open('converted_mat/'+filename+'.txt', 'w') as writefile:
			writefile.write(str(data))
		# pprint (data)
