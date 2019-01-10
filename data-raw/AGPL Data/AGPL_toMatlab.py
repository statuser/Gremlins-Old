import pandas as pd

data = pd.read_csv("AGPL-Training.csv",index_col = 0)
data['V13'] = data['V13']/1000

print data.head()
data.to_csv('AgplData.txt',header=None,index=None,sep='\t')

