import numpy as np
import pandas as pd
import warnings
import re
from sklearn import pipeline
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn.compose import ColumnTransformer
warnings.filterwarnings('ignore')

GPL = pd.read_csv("GPL201-30390.csv")
GSE1739 = pd.read_csv("GSE1739.csv")

GPL.info()
GSE1739.info()

GSE1739.insert(15, 'Gene Symbol', GPL['Gene Symbol'])
for key, values in GSE1739.iteritems():
    print(key, values)
GSE1739 = GSE1739.dropna(axis='index', subset='Gene Symbol')
tmp = GSE1739['Gene Symbol']

regex = r"\d{2}(:)\d{2}(:)\d{2}"
def CheckStr(s):
    if re.findall(regex, s):
        return 1
    else:
        return 0
    
matched_index = []
for key, value in tmp.items():
    if CheckStr(value) == 1:
        matched_index.append(key)
        
print(matched_index)
matched_index = np.array(matched_index)
print("Date rows: ", matched_index)
GSE1739.loc[matched_index,:]
GSE1739 = GSE1739.drop(axis='index', index=matched_index)

tmp = GSE1739['Gene Symbol'].values.tolist()
sym = []
for s in tmp:
    str = '' ''
    s = s.split('///')
    s = str.join(s)
    sym.append(s)

S = {'Gene Symbol':sym}
GS = pd.DataFrame(S)

GSE1739 = GSE1739.drop(columns='Gene Symbol')
GSE1739 = pd.concat([GSE1739, GS], axis=1, join='inner')

print(GSE1739.isna().sum())
GSE1739.info()

num_cols = list(GSE1739.select_dtypes(include = 'float64').columns)
norm_pipeline = Pipeline([('scaler', MinMaxScaler())])
preprocess_pipeline = ColumnTransformer([('numcols', norm_pipeline, num_cols)])
tmp = preprocess_pipeline.fit_transform(GSE1739)

gid = np.array(GSE1739['ID'].tolist())
gs = np.array(GSE1739['Gene Symbol'].tolist())
print(gid, gs)

cols = ['GSM30361', 'GSM30362', 'GSM30363', 'GSM30364',
        'GSM30365', 'GSM30366', 'GSM30367', 'GSM30368',
        'GSM30369', 'GSM30370', 'GSM30371', 'GSM30372',
        'GSM30373', 'GSM30374', 'ID', 'Gene Symbol']
tmp = np.column_stack((gid, tmp, gs))
tmp = pd.DataFrame(tmp, columns=cols)

GSE1739 = tmp
print(GSE1739)

# GSE1739.to_csv('GSE1739(Processed).csv')