#!/usr/bin/env python
# coding: utf-8

# # Input (VCF+SEI)

# In[1]:


import numpy as np
import pandas as pd
combined4colbasedf = pd.read_csv("/data/projects/nanopore/RepeatExpansion/TR_downstreamAnalysis/vcf/combined4colbase.vcf",sep="\t", header=None, names=['chrX','1','2','3','4'])
combined4colbasedf
df=pd.read_csv("/data/projects/nanopore/RepeatExpansion/TR_downstreamAnalysis/tsv/combined4colbase.ref_combined.tsv", sep="\t",low_memory=False)
df = df.drop_duplicates(subset=["1"])
print(df)

frames= [df, combined4colbasedf]
finalinput =pd.merge(right=combined4colbasedf, left=df, on=["1","2"])
dfinput= finalinput
display(dfinput)


dfinput.columns = [ "0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","chromosome","window","basepair"]
dfinput


display(dfinput)
dfinput = dfinput.join(dfinput['window'].str.split('_', expand=True).rename (columns={0:'TR_id', 1:'Win_num'}))
dfinput['Win_num']=dfinput['Win_num'].astype(int)
dfinput['TR_id']=dfinput['TR_id'].astype(int)
subset= dfinput[dfinput.Win_num<=50] 
subset= subset.sort_values(by=['TR_id','Win_num'])
subset

dfUpstreamdropwinnum = subset.drop_duplicates(subset=["TR_id","Win_num"], keep="first") 
print(dfUpstreamdropwinnum)

dfDownstreamdropwinnum = subset.drop_duplicates(subset=["TR_id","Win_num"], keep="last") 
print(dfDownstreamdropwinnum)


# # Adding the Column to iterate thourgh, for Matrix Flattening

# In[2]:


dfDownstreamdropwinnum=(dfDownstreamdropwinnum.reset_index(drop=True))
dfDownwinnum = dfDownstreamdropwinnum.drop(columns=['0', '1','2','3','4','5','6','7','8', 'chromosome','window','basepair'])


result = []
i = 0
for j in range(len(dfDownwinnum["TR_id"])):
   
    
    if j == len(dfDownwinnum["TR_id"])-1:
        result.append(i)
        
    elif dfDownwinnum["TR_id"].iloc[j-1] != dfDownwinnum["TR_id"].iloc[j]:
        result.append(i+1)
        i=i+1
          # if j ==0 append (i) 
    else:
        result.append(i)

dfDownwinnum["Result"] = result  
print(dfDownwinnum)


# In[3]:


dfUpstreamdropwinnum=(dfUpstreamdropwinnum.reset_index(drop=True))
dfUpwinnum = dfUpstreamdropwinnum.drop(columns=['0', '1','2','3','4','5','6','7','8', 'chromosome','window','basepair'])

result = []
i = 0
for j in range(len(dfUpwinnum["TR_id"])):
   
    
    if j == len(dfUpwinnum["TR_id"])-1:
        result.append(i)
        
    elif dfUpwinnum["TR_id"].iloc[j-1] != dfUpwinnum["TR_id"].iloc[j]:
        result.append(i+1)
        i=i+1
          # if j ==0 append (i) 
    else:
        result.append(i)

dfUpwinnum["Result"] = result  
print(dfUpwinnum)


# # Matrix Flattening 

# In[4]:




def condense_df(df, tr_id):
    #UpstreamMat = df.loc[df["TR_id"] == str(tr_id)].copy()
    DownstreamMatwinum = df.loc[df["Result"] == (tr_id)].copy()
    DownstreamMatwinum.drop("Result", axis=1,inplace=True)
    DownstreamMatwinum.drop("TR_id", axis=1,inplace=True)
    DownstreamMatwinum.drop("Win_num", axis=1,inplace=True)
    arrDown = DownstreamMatwinum.to_numpy().flatten(order='F')
    return arrDown


colsDOWN = list()
colsDOWN = dfDownwinnum.columns.tolist()
#cols_newUP = [colsUP[-1]]
print(type(colsDOWN))
print(len(colsDOWN))
#cols_newUP.extend(colsUP[0:]) 
DownstreamMatrix=dfDownwinnum
#print(UpstreamMatrix.columns)

DownstreamMatwinum = np.zeros(shape=(len(set(DownstreamMatrix["Result"].tolist())), len(condense_df(DownstreamMatrix, 1))))

print(DownstreamMatwinum.shape)
cnt = 0
failed_ids= []

for i in list(set(DownstreamMatrix["Result"].tolist())):
    cnt +=1
    try:
        DownstreamMatwinum [int (i)-1,:] = condense_df(DownstreamMatrix, i)
    except:
        failed_ids.append(i) 



DownstreamMatwinum_copy = DownstreamMatwinum
failed_ids = [i-1 for i in list(map(int,failed_ids))]
DownstreamMatdropwinum = np.delete (DownstreamMatwinum_copy, failed_ids, axis = 0)
np.save("DownstreamMatwinum",DownstreamMatdropwinum)
print(DownstreamMatdropwinum)


# In[5]:




def condense_df(df, tr_id):
    #UpstreamMat = df.loc[df["TR_id"] == str(tr_id)].copy()
    UpstreamMatwinum = df.loc[df["Result"] == (tr_id)].copy()
    UpstreamMatwinum.drop("Result", axis=1,inplace=True)
    UpstreamMatwinum.drop("TR_id", axis=1,inplace=True)
    UpstreamMatwinum.drop("Win_num", axis=1,inplace=True)
    arrUp = UpstreamMatwinum.to_numpy().flatten(order='F')
    return arrUp


colsUp = list()
colsUp = dfUpwinnum.columns.tolist()
#cols_newUP = [colsUP[-1]]
print(type(colsUp))
print(len(colsUp))
#cols_newUP.extend(colsUP[0:]) 
UpstreamMatrix=dfUpwinnum
#print(UpstreamMatrix.columns)

UpstreamMatwinum = np.zeros(shape=(len(set(UpstreamMatrix["Result"].tolist())), len(condense_df(UpstreamMatrix, 1))))

print(UpstreamMatwinum.shape)
cnt = 0
failed_ids= []

for i in list(set(UpstreamMatrix["Result"].tolist())):
    cnt +=1
    try:
        UpstreamMatwinum [int (i)-1,:] = condense_df(UpstreamMatrix, i)
    except:
        failed_ids.append(i) 



UpstreamMatwinum_copy = UpstreamMatwinum
failed_ids = [i-1 for i in list(map(int,failed_ids))]
UpstreamMatdropwinum = np.delete (UpstreamMatwinum_copy, failed_ids, axis = 0)
np.save("UpstreamMatwinum",UpstreamMatdropwinum)
print(UpstreamMatdropwinum)


# In[6]:


log10downstreamdropwin = np.log10(DownstreamMatdropwinum)
print(log10downstreamdropwin)


# In[7]:


log10upstreamdropwin = np.log10(UpstreamMatdropwinum)
print(log10upstreamdropwin)


# # DOWNSTREAM PCA

# In[8]:


from sklearn import datasets 
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

Downstreamdropwin= pd.DataFrame(log10downstreamdropwin)
print(Downstreamdropwin)

Downstreamdropwin.to_csv("Downstreamdropwin", index= None)
Downstreamdropwin = pd.read_csv("Downstreamdropwin")

Downstreamdropwin
x = Downstreamdropwin.values
y = Downstreamdropwin.values
x = StandardScaler().fit_transform(x)

from sklearn.decomposition import PCA
pcadfdownwinnum = PCA(n_components=8)
pcadfdownwinnummatrix = pcadfdownwinnum.fit_transform(x)

plt.bar(x=range(8), height= pcadfdownwinnum.explained_variance_ratio_)

plt.show()

sum(pcadfdownwinnum.explained_variance_ratio_)


principaldownstreamwinnum = pd.DataFrame (data = pcadfdownwinnummatrix, columns = ['a', 'b','c','d','e','f','g','h'])
plt.scatter(principaldownstreamwinnum['a'], principaldownstreamwinnum['b'], c='green')
plt.show()


# # UPSTREAM PCA

# In[9]:




Upstreamdropwin= pd.DataFrame(log10upstreamdropwin)
print(Upstreamdropwin)

Upstreamdropwin.to_csv("Upstreamdropwin", index= None)
Upstreamdropwin = pd.read_csv("Upstreamdropwin")

Upstreamdropwin
x = Upstreamdropwin.values
y = Upstreamdropwin.values
x = StandardScaler().fit_transform(x)

from sklearn.decomposition import PCA
pcadfupwinnum = PCA(n_components=8)
pcadfupwinnummatrix = pcadfupwinnum.fit_transform(x)

plt.bar(x=range(8), height= pcadfupwinnum.explained_variance_ratio_)

plt.show()

sum(pcadfupwinnum.explained_variance_ratio_)

principalupstreamwinnum = pd.DataFrame (data = pcadfupwinnummatrix, columns = ['a', 'b','c','d','e','f','g','h'])
plt.scatter(principalupstreamwinnum['a'], principalupstreamwinnum['b'], c='green')
plt.show()


# In[15]:


xup= StandardScaler().fit_transform(x)


# In[20]:


Downstreamdropwin


# # TSNE DownStream 

# In[11]:


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas

import sklearn

from sklearn.manifold import TSNE
tsne = TSNE(n_components=2,
              random_state=12)


n_components=2
tsne = TSNE(n_components=n_components,
              perplexity=1000,
              random_state=12)
Z= Downstreamdropwin.values
X_2d1 = tsne.fit_transform(Z)



X_2d1 = pd.DataFrame (data = X_2d1, columns = ['a', 'b'])
plt.scatter(X_2d1 ['a'], X_2d1 ['b'], c='green')
plt.title('With perplexity = 1000, tsne for downstream Matrix after dropping winnum')
plt.show()


n_components=2
tsne = TSNE(n_components=n_components,
              perplexity=50,
              random_state=12)
Z= Downstreamdropwin.values
X_2d2 = tsne.fit_transform(Z)



X_2d2 = pd.DataFrame (data = X_2d2, columns = ['a', 'b'])
plt.scatter(X_2d2 ['a'], X_2d2 ['b'], c='green')
plt.title('With perplexity = 50, tsne for downstream Matrix after dropping winnum')
plt.show()



n_components=2
tsne = TSNE(n_components=n_components,
              perplexity=5,
              random_state=12)
Z= Downstreamdropwin.values
X_2d = tsne.fit_transform(Z)



X_2d = pd.DataFrame (data = X_2d, columns = ['a', 'b'])
plt.scatter(X_2d ['a'], X_2d ['b'], c='green')
plt.title('With perplexity = 50, tsne for downstream Matrix after dropping winnum')
plt.show()


# In[ ]:


X_2d 


# # TSNE Upstream

# In[12]:


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas

import sklearn

from sklearn.manifold import TSNE
tsne = TSNE(n_components=2,
              random_state=12)


n_components=2
tsne = TSNE(n_components=n_components,
              perplexity=1000,
              random_state=12)
Z= Upstreamdropwin.values
X_2u = tsne.fit_transform(Z)



X_2u = pd.DataFrame (data = X_2u, columns = ['a', 'b'])
plt.scatter(X_2u ['a'], X_2u ['b'], c='green')
plt.title('With perplexity = 1000, tsne for upstream Matrix after dropping winnum')
plt.show()


n_components=2
tsne = TSNE(n_components=n_components,
              perplexity=50,
              random_state=12)
Z= Upstreamdropwin.values
X_2u1 = tsne.fit_transform(Z)



X_2u1 = pd.DataFrame (data = X_2u1, columns = ['a', 'b'])
plt.scatter(X_2u1 ['a'], X_2u1 ['b'], c='green')

plt.title('With perplexity = 50, tsne for Upstreamm Matrix after dropping winnum')
plt.show()

n_components=2
tsne = TSNE(n_components=n_components,
              perplexity=5,
              random_state=12)
Z= Upstreamdropwin.values
X_2u2 = tsne.fit_transform(Z)



X_2u2 = pd.DataFrame (data = X_2u2, columns = ['a', 'b'])
plt.scatter(X_2u2 ['a'], X_2u2 ['b'], c='green')
plt.title('With perplexity = 5, tsne for Upstream Matrix after dropping winnum')
plt.show()


# In[ ]:


X_2u2


# # Upstream coloring (Not usable but keeping it as a reference)

# In[ ]:


principalupstreamwinnum[['a','b','c','e','f','g','h']] 


## hue corresponds to the value
## write the function that would do everything 
## dark - most similar features 
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principalupstreamwinnum, x="a", y="b", hue="b")


# # Downstream coloring (Not Quite Usable, but Keeping it as a reference)

# In[ ]:


principaldownstreamwinnum[['a','b','c','e','f','g','h']] 


## hue corresponds to the value
## write the function that would do everything 
## dark - most similar features 
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue="c")


# # DF Winnmum UP

# In[13]:


def sum_df(df, tr_id, seqclass):
    #print (df)
    dfUpwinnum = df.loc[df["Result"] == (tr_id)].copy()
    sumTRUP =dfUpwinnum[str(9+seqclass)].sum()
    
    return sumTRUP



# In[14]:


colsUP = list()
colsUP = dfUpwinnum.columns.tolist()
print(type(colsUP))
print(len(colsUP))
#cols_newDOWN.extend(colsDOWN[0:]) 

#print(DownstreamMatrix.columns)


cnt = 0
sumTRUP= {'PC1': [], 'E3': [], 'E4': [], 'HET1': [], 'E8': [], 'HET2': [], 'E9': [], 'HET3':[], 'PC4' : [], 'P': [], 'CTCF' : [], 'E10' : [], 'HET4': []}
translater = {'PC1': 0, 'E3': 7, 'E4': 9, 'HET1': 11, 'E8': 17, 'HET2': 23,'E9': 26, 'HET3':29, 'PC4' :34, 'P': 25, 'CTCF' : 27, 'E10' : 30, 'HET4': 32}

for i in list(set(dfUpwinnum["Result"].tolist())):
    for key in translater.keys():
        sumTRUP[key].append(sum_df(dfUpwinnum,i, translater[key]))

dfUpwinnumdrop = dfUpwinnum


# In[16]:


for key in translater.keys():
    plt.hist(sumTRUP[key])
    plt.title(key)
    plt.show()


# In[17]:


for key in translater.keys():
    plt.plot(sumTRUP[key])
    plt.title(key)
    plt.show()


# In[18]:


for key in translater.keys():

    X_2u1[key] = sumTRUP[key]
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.scatterplot(data=X_2u1, x="a", y="b", hue=key)
    plt.title(key)
    plt.show()


# In[19]:


for key in translater.keys():

    X_2u2[key] = sumTRUP[key]
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.scatterplot(data=X_2u2, x="a", y="b", hue=key)
    plt.title(key)
    plt.show()
    
    


# # UPStream PCA coloring by SumTR

# In[30]:


import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


## if you have any existing df and want to add columns (the same # of rows,take the data) = add as a column 
## 3 important, take the data 3 columns ( easy way to create the columns) 
principalupstreamwinnum[key] = sumTRUP[key]

## hue corresponds to the value
## write the function that would do everything 
## dark - most similar features 
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principalupstreamwinnum, x="a", y="b", hue=key)


# In[ ]:


dfUpwinnum


# In[ ]:


dfUpwinnumdrop


# In[ ]:


fig,ax = plt.subplots(-(-len(translater.keys())//3),3)

for i,key in enumerate(translater.keys()):   
    
    ax[i%4][i%2].hist(sumTR[key])
    ax[i%4][i%2].set_title(key)
    

    #yaxis: num TR xaxis: sumval
    #hist: how common are the values 


# In[27]:


for key in translater.keys():

    principalupstreamwinnum[key] = sumTRUP[key]
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.scatterplot(data=principalupstreamwinnum, x="a", y="b", hue=key)
    plt.title(key)
    plt.show()

    


# # Subset by more clustered region Tsne

# In[ ]:


dfUpwinnum[dfUpwinnum.Result.isin(X_2d[X_2d.P > 750].index+1)]


# In[ ]:


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas

import sklearn

from sklearn.manifold import TSNE
tsne = TSNE(n_components=2,
              random_state=12)


n_components=2
tsne = TSNE(n_components=n_components,
              perplexity=1000,
              random_state=12)
Z= dfUpwinnum.values
X_2d = tsne.fit_transform(Z)



X_2d = pd.DataFrame (data = X_2d, columns = ['a', 'b'])
plt.scatter(X_2d ['a'], X_2d ['b'], c='green')
plt.title('With perplexity = 1000, tsne for dfUpwinnum')
plt.show()


n_components=2
tsne = TSNE(n_components=n_components,
              perplexity=50,
              random_state=12)
Z= dfUpwinnum.values
X_2d = tsne.fit_transform(Z)



X_2d = pd.DataFrame (data = X_2d, columns = ['a', 'b'])
plt.scatter(X_2d ['a'], X_2d ['b'], c='green')
plt.title('With perplexity = 50, tsne for dfUpwinnum')
plt.show()


# # DF Winnmum DOWN

# In[15]:


def sum_df(df, tr_id, seqclass):
    #print (df)
    dfDownwinnum = df.loc[df["Result"] == (tr_id)].copy()
    sumTR =dfDownwinnum[str(9+seqclass)].sum()
    
    return sumTR



colsDown = list()
colsDown = dfDownwinnum.columns.tolist()
print(type(colsDown))
print(len(colsDown))
#cols_newDOWN.extend(colsDOWN[0:]) 

#print(DownstreamMatrix.columns)


cnt = 0
sumTR= {'PC1': [], 'E3': [], 'E4': [], 'HET1': [], 'E8': [], 'HET2': [], 'E9': [], 'HET3':[], 'PC4' : [], 'P': [], 'CTCF' : [], 'E10' : [], 'HET4': []}
translater = {'PC1': 0, 'E3': 7, 'E4': 9, 'HET1': 11, 'E8': 17, 'HET2': 23,'E9': 26, 'HET3':29, 'PC4' :34, 'P': 25, 'CTCF' : 27, 'E10' : 30, 'HET4': 32}

for i in list(set(dfDownwinnum["Result"].tolist())):
    for key in translater.keys():
        sumTR[key].append(sum_df(dfDownwinnum,i, translater[key]))

dfDownwinnumdrop = dfDownwinnum


# In[32]:


for key in translater.keys():
    plt.plot(sumTR[key])
    plt.title(key)
    plt.show()


# In[33]:


for key in translater.keys():

    X_2d[key] = sumTR[key]
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.scatterplot(data=X_2d, x="a", y="b", hue=key)
    plt.title(key)
    plt.show()


# # DOWNStream PCA coloring by SumTR

# In[34]:


for key in translater.keys():

    principaldownstreamwinnum[key] = sumTR[key]
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue=key)
    plt.title(key)
    plt.show()
    


# In[35]:


for key in translater.keys():

    X_2d[key] = sumTR[key]
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.scatterplot(data=X_2d, x="a", y="b", hue=key)
    plt.title(key)
    plt.show()


# In[36]:


for key in translater.keys():

    principaldownstreamwinnum[key] = sumTR[key]
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue=key)
    plt.title(key)
    plt.show()


# # most Contirbuting 3 (UPstream)

# In[ ]:


#pca.component , upstream first 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


pca = PCA()
df_pca = pca.fit_transform(Upstreamdropwin)

most_important3_differentmethod=np.abs(pca.components_).argsort()[::-1][:3]

mostimportant3= most_important3_differentmethod[:,0]
mostimportant3= list(mostimportant3)


# In[ ]:


most_important3_differentmethod


# In[ ]:


Upstreamdropwin.columns = Upstreamdropwin.columns.astype(int) 
Upstreamdropwin.columns 


# In[ ]:


mostimportantsubset=Upstreamdropwin[mostimportant3]
mostimportantsubset


# In[ ]:


from sklearn.decomposition import PCA
pcadfupstreamMatreal = PCA(n_components=2)
principalComponentsdfupstreamMatreal = pcadfupstreamMatreal.fit_transform(Upstreamdropwin)

plt.bar(x=range(2), height= pcadfupstreamMatreal.explained_variance_ratio_)

plt.show()

sum(pcadfupstreamMatreal.explained_variance_ratio_)

principalupstreamDfreal = pd.DataFrame (data = principalComponentsdfupstreamMatreal, columns = ['a', 'b'])
plt.scatter(principalupstreamDfreal['a'], principalupstreamDfreal['b'], c='purple')
plt.show()


# In[ ]:


import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


## if you have any existing df and want to add columns (the same # of rows,take the data) = add as a column 
## 3 important, take the data 3 columns ( easy way to create the columns) 
principalupstreamwinnum[['mostimportant_1','most_important_2','mostimportant_3']] = mostimportantsubset

## hue corresponds to the value
## write the function that would do everything 
## dark - most similar features 
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principalupstreamwinnum, x="a", y="b", hue="mostimportant_1")


# # most Contirbuting 3 (DOWNstream)

# In[ ]:


#pca.component , upstream first 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


pca = PCA()
df_pca = pca.fit_transform(Downstreamdropwin)

most_important3ones_differentmethod=np.abs(pca.components_).argsort()[::-1][:3]

mostimportant3ones= most_important3ones_differentmethod[:,0]
mostimportant3ones= list(mostimportant3)


# In[ ]:


Downstreamdropwin.columns = Downstreamdropwin.columns.astype(int) 
Downstreamdropwin.columns 


# In[ ]:


mostimportantsubsetDown=Downstreamdropwin[mostimportant3ones]
mostimportantsubsetDown


# In[ ]:


from sklearn.decomposition import PCA
pcadfdownstreamMatreal = PCA(n_components=2)
principalComponentsdfdownstreamMatreal = pcadfdownstreamMatreal.fit_transform(Downstreamdropwin)

plt.bar(x=range(2), height= pcadfdownstreamMatreal.explained_variance_ratio_)

plt.show()

sum(pcadfdownstreamMatreal.explained_variance_ratio_)

principaldownstreamDfreal = pd.DataFrame (data = principalComponentsdfdownstreamMatreal, columns = ['a', 'b'])
plt.scatter(principaldownstreamDfreal['a'], principaldownstreamDfreal['b'], c='green')
plt.show()


# In[ ]:


import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


## if you have any existing df and want to add columns (the same # of rows,take the data) = add as a column 
## 3 important, take the data 3 columns ( easy way to create the columns) 
principaldownstreamwinnum[['mostimportant_1','most_important_2','mostimportant_3']] = mostimportantsubsetDown

## hue corresponds to the value
## write the function that would do everything 
## dark - most similar features 
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue="mostimportant_1")


# # DF winnum subset by SumTR data  (according to P) /Upstream

# In[ ]:


SubdfUpwim=dfUpwinnum[dfUpwinnum.Result.isin(X_2d[X_2d.P > 750].index+1)]
SubdfUpwim


# In[ ]:


newSubdfUpwim=(SubdfUpwim.reset_index(drop=True))


result = []
i = 0
for j in range(len(newSubdfUpwim["TR_id"])):
   
    
    if j == len(newSubdfUpwim["TR_id"])-1:
        result.append(i)
        
    elif newSubdfUpwim["TR_id"].iloc[j-1] != newSubdfUpwim["TR_id"].iloc[j]:
        result.append(i+1)
        i=i+1
          # if j ==0 append (i) 
    else:
        result.append(i)

newSubdfUpwim["Seq"] = result  
print(newSubdfUpwim)


# In[ ]:


newSubdfUpwim.to_csv('newSubdfUpwim.tsv',sep='\t',index = False)

def condense_df(df, tr_id):
    #UpstreamMat = df.loc[df["TR_id"] == str(tr_id)].copy()
    newSubdfUpwinums= df.loc[df["Seq"] == (tr_id)].copy()
    newSubdfUpwinums.drop("Result", axis=1,inplace=True)
    newSubdfUpwinums.drop("TR_id", axis=1,inplace=True)
    newSubdfUpwinums.drop("Win_num", axis=1,inplace=True)
    newSubdfUpwinums.drop("Seq", axis=1,inplace=True)
    arrNew = newSubdfUpwinums.to_numpy().flatten(order='F')
    return arrNew


colsUPs = list()
colsUPs = newSubdfUpwim.columns.tolist()
#cols_newUP = [colsUP[-1]]
print(type(colsUPs))
print(len(colsUPs))
#cols_newUP.extend(colsUP[0:]) 
newSubdfUpwims=newSubdfUpwim
print(newSubdfUpwims.shape)
#print(UpstreamMatrix.columns)

newSubdfUpwinums = np.zeros(shape=(len(set(newSubdfUpwims["Seq"].tolist())), len(condense_df(newSubdfUpwims, 1))))

print(newSubdfUpwinums.shape)
cnt = 0
failed_ids= []

for i in list(set(newSubdfUpwims["Seq"].tolist())):
    cnt +=1
    try:
        newSubdfUpwinums[int (i)-1,:] = condense_df(newSubdfUpwims, i)
    except:
        failed_ids.append(i) 



newSubdfUpwims_copy = newSubdfUpwinums
failed_ids = [i-1 for i in list(map(int,failed_ids))]
newSubdfUpwinumms = np.delete (newSubdfUpwims_copy, failed_ids, axis = 0)
np.save("newSubdfUpwinums",newSubdfUpwinumms)
print(newSubdfUpwinumms)


# In[ ]:


log10newSubdfUpwinumms = np.log10(newSubdfUpwinumms)
print(log10newSubdfUpwinumms)


# In[ ]:


from sklearn import datasets 
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

log10newSubdfUpwinumms= pd.DataFrame(log10newSubdfUpwinumms)
print(log10newSubdfUpwinumms)

log10newSubdfUpwinumms.to_csv("log10newSubdfUpwinumms", index= None)
log10newSubdfUpwinumms = pd.read_csv("log10newSubdfUpwinumms")

log10newSubdfUpwinumms
x = log10newSubdfUpwinumms.values
y = log10newSubdfUpwinumms.values
x = StandardScaler().fit_transform(x)

from sklearn.decomposition import PCA
pcadflog10newSubdfUpwinumms = PCA(n_components=8)
pcadflog10newSubdfUpwinummsmatrix = pcadflog10newSubdfUpwinumms.fit_transform(x)

plt.bar(x=range(8), height= pcadflog10newSubdfUpwinumms.explained_variance_ratio_)

plt.show()

sum(pcadflog10newSubdfUpwinumms.explained_variance_ratio_)

principalpcadflog10newSubdfUpwinumms = pd.DataFrame (data = pcadflog10newSubdfUpwinummsmatrix, columns = ['a', 'b','c','d','e','f','g','h'])
plt.scatter(principalpcadflog10newSubdfUpwinumms['a'], principalpcadflog10newSubdfUpwinumms['b'], c='green')
plt.show()


# In[65]:


for key in translater.keys():

    principalpcadflog10newSubdfUpwinumms[key] = sumTR[key]
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.scatterplot(data=principalpcadflog10newSubdfUpwinumms, x="a", y="b", hue=key)
    plt.title(key)
    plt.show()
    


# # Nearest Exon Coloring

# In[ ]:


# input 


# In[ ]:


nearExonTRID = pd.read_csv("sorteddropdup.bed",sep="\t", header=None,names=['1','2','TR_id'])
nearExonTRID

nearExon = pd.read_csv("nearestExon.bed",sep="\t", header=None,names=['chrX','1','2','3','4','5','6'])
nearExon

dropExon = nearExon.drop_duplicates(subset=["2","3"],keep="first")
dropExon


frames= [dropExon, nearExonTRID]
finalExoninput =pd.merge(right=nearExonTRID, left=dropExon, on=["1","2"])

display(finalExoninput)


# In[ ]:


np.where(finalExoninput['6']>5000, True, False)


# # UPstream Exon PCA/TSNE data

# In[ ]:


plt.hist(finalExoninput['6'])


# In[ ]:


principalupstreamwinnum['exon'] = np.where((finalExoninput['6']>-5000)&(finalExoninput['6']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principalupstreamwinnum, x="a", y="b", hue='exon')
plt.title(key)
plt.show()


# # DOWNStream Exon PCA/TSNE data

# In[ ]:


principaldownstreamwinnum['exon'] = np.where((finalExoninput['6']>-5000)&(finalExoninput['6']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue='exon')
plt.title(key)
plt.show()


# In[ ]:


5


# In[ ]:


finalExoninput['6'].value_counts()


# # transcription start site Coloring

# In[ ]:


TSS_coord = pd.read_csv("/data/projects/nanopore/RepeatExpansion/coordinates/TSS_coords.bed",sep="\t", header=None, names=['chrX','0','1','2','3'])
TSS_coord


# In[ ]:


import numpy as np
import pandas as pd
TSS_merged = pd.read_csv("/data/projects/nanopore/RepeatExpansion/TR_downstreamAnalysis/TSS/testTSSmerge.bed", sep="\t", header=None, names=['chrX','1','2','chrX1','4','5','6'])
TSS_merged 


# In[ ]:


dropTSS_merged=TSS_merged.drop_duplicates(subset=["1","2"],keep="first") 
dropTSS_merged


# # Upstream PCA/TSNE with TSS

# In[ ]:


X_2u['TSS'] = np.where((dropTSS_merged['6']>-5000)&(dropTSS_merged['6']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u, x="a", y="b", hue='TSS')
plt.title('Upstream TSNE with TSS')
plt.show()


# In[ ]:


X_2u1['TSS'] = np.where((dropTSS_merged['6']>-5000)&(dropTSS_merged['6']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u1, x="a", y="b", hue='TSS')
plt.title('Upstream TSNE with TSS')
plt.show()


# In[ ]:


X_2u2['TSS'] = np.where((dropTSS_merged['6']>-5000)&(dropTSS_merged['6']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u2, x="a", y="b", hue='TSS')
plt.title('Upstream TSNE with TSS')
plt.show()


# In[ ]:


principalupstreamwinnum['TSS'] = np.where((dropTSS_merged['6']>-5000)&(dropTSS_merged['6']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principalupstreamwinnum, x="a", y="b", hue='TSS')
plt.title("TSS upstram")
plt.show()


# # Downstream PCA/TSNE with TSS data

# In[ ]:


principaldownstreamwinnum['TSS'] = np.where((dropTSS_merged['6']>-5000)&(dropTSS_merged['6']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue='TSS')
plt.title("PCA and TSS DOWNSTREAM")
plt.show()


# In[ ]:


X_2d['TSS'] = np.where((dropTSS_merged['6']>-5000)&(dropTSS_merged['6']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='TSS')
plt.title("TSNE and TSS downstream")
plt.show()


# # CTCF Coloring

# In[ ]:


import pandas as pd
import numpy as np


# In[ ]:


CTCF = pd.read_csv("/data/projects/nanopore/RepeatExpansion/TR_downstreamAnalysis/CTCF/droppedchrXCTCFClosest.bed", sep="\t", header=None, names=['1','2','chrx','3','4','CTCF','closet'])
CTCF


# In[ ]:


CTCF_merged=CTCF.drop_duplicates(subset=["1","2"],keep="first") 
CTCF_merged


# In[ ]:


X_2u['CTCF'] = np.where((CTCF_merged['closet']>-5000)&(CTCF_merged['closet']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u, x="a", y="b", hue='CTCF')
plt.title('Upstream TSNE with CTCF')
plt.show()


# In[ ]:


X_2d['CTCF'] = np.where((CTCF_merged['closet']>-5000)&(CTCF_merged['closet']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='CTCF')
plt.title('Downstream TSNE with CTCF')
plt.show()


# In[ ]:


X_2u2['CTCF'] = np.where((CTCF_merged['closet']>-5000)&(CTCF_merged['closet']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u2, x="a", y="b", hue='CTCF')
plt.title('Upstream TSNE with CTCF')
plt.show()


# 

# # region length Coloring

# In[ ]:


regionlen = pd.read_csv("/data/projects/nanopore/RepeatExpansion/TR_downstreamAnalysis/closetregionlen.bed", sep="\t", header=None, names=['chrx','1','2','chrx1','3','4','length','tr_id','null'])
regionlen


# In[ ]:


X_2d['length'] =regionlen['length']
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='length')
plt.title("tsne for downstream, TSNE")
plt.show()
    


# In[ ]:


X_2u['length'] =regionlen['length']
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u, x="a", y="b", hue='length')
plt.title("upstream tsne with Region Length")
plt.show()


# In[ ]:


X_2u2['length'] =regionlen['length']
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u2, x="a", y="b", hue='length')
plt.title("upstream tsne with Region Length")
plt.show()


# In[ ]:


principaldownstreamwinnum['length'] =regionlen['length']
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue='length')
plt.title("downstream pca with Region Length")
plt.show()


# In[ ]:


principalupstreamwinnum['length'] =regionlen['length']
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principalupstreamwinnum, x="a", y="b", hue='length')
plt.title("upstream pca with Region length")
plt.show()


# # Repeat Masker/the repeatmasker classes corresponde to all transposon and retrotransposon classes

# In[ ]:


import numpy as np
import pandas as pd
repeatmasker = pd.read_csv("/data/projects/nanopore/RepeatExpansion/TR_downstreamAnalysis/closestsortrepeatchrX.bed", sep="\t", header=None, names=['chrX','1','2','chrX1','3','4','class','distance'])
repeatmasker


# In[ ]:


droprepeatmasker=repeatmasker.drop_duplicates(subset=["1","2"],keep="first") 
droprepeatmasker


# In[ ]:


principalupstreamwinnum['distance'] =droprepeatmasker['distance']
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principalupstreamwinnum, x="a", y="b", hue='distance')
plt.title("upstream pca with repeatmasker by distance")
plt.show()


# # filter by distance (Repeat Masker)

# In[ ]:


X_2u['distance'] = np.where((droprepeatmasker['distance']>-200)&(droprepeatmasker['distance']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u, x="a", y="b", hue='distance')
plt.title('Upstream TSNE with Repeat Masker Distance')
plt.show()


# In[ ]:


X_2d['distance'] = np.where((droprepeatmasker['distance']<200)&(droprepeatmasker['distance']>0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='distance')
plt.title('Downstream TSNE with Repeat Masker Distance')
plt.show()


# In[ ]:


X_2d['distance'] = np.where((droprepeatmasker['distance']>-200)&(droprepeatmasker['distance']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='distance')
plt.title('Downstream TSNE with Repeat Masker Distance')
plt.show()


# In[ ]:


principalupstreamwinnum['distance'] = np.where((droprepeatmasker['distance']>-200)&(droprepeatmasker['distance']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principalupstreamwinnum, x="a", y="b", hue='distance')
plt.title('UPstream pca with Repeat Masker Distance')
plt.show()


# In[ ]:


principaldownstreamwinnum['distance'] = np.where((droprepeatmasker['distance']>-200)&(droprepeatmasker['distance']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue='distance')
plt.title('Downstream pca with Repeat Masker Distance')
plt.show()


# In[ ]:


principaldownstreamwinnum['distance'] = np.where((droprepeatmasker['distance']<200)&(droprepeatmasker['distance']>0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue='distance')
plt.title('Downstream pca with Repeat Masker Distance')
plt.show()


# # assign to Classes such as LINE, Alu, HERV, MER, LTR, L2 (Repeat Masker)

# In[ ]:


principalupstreamwinnum['class'] =droprepeatmasker['class']
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principalupstreamwinnum, x="a", y="b", hue='class')
plt.title("upstream pca with repeatmasker")
plt.show()


# In[ ]:


# pca downstream 
principalupstreamwinnum['class'] =droprepeatmasker['class']
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principalupstreamwinnum, x="a", y="b", hue='class')
plt.title("upstream pca with repeatmasker")
plt.show()

# pca upstream
principaldownstreamwinnum['class'] =droprepeatmasker['class']
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue='class')
plt.title("downstream pca with Repeatmasker")
plt.show()

#tsne down
X_2d['class']=np.where(droprepeatmasker['class'].str.contains("LTR"),"LTR",droprepeatmasker['class'])
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='class')
plt.title('Downstream TSNE with repeatmasker')
plt.show()

#tsne up
X_2u['class']=np.where(droprepeatmasker['class'].str.contains("LTR"),"LTR","MER")
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='class')
plt.title('Upstream TSNE with repeatmasker')
plt.show()

X_2u1['class']=np.where(droprepeatmasker['class'].str.contains("LTR"),"LTR","MER")
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u1, x="a", y="b", hue='class')
plt.title('Upstream TSNE with repeatmasker')
plt.show()

X_2u2['class']=np.where(droprepeatmasker['class'].str.contains("LTR"),"LTR","MER")
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u2, x="a", y="b", hue='class')
plt.title('Upstream TSNE with repeatmasker')
plt.show()


# In[ ]:


X_2d['class']=np.where(droprepeatmasker['class'].str.contains("LTR"),"LTR",droprepeatmasker['class'])
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='class')
plt.title('Downstream TSNE with repeatmasker')
plt.show()


# In[ ]:


X_2u['class']=np.where(droprepeatmasker['class'].str.contains("LTR"),"LTR","MER")
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='class')
plt.title('Upstream TSNE with repeatmasker')
plt.show()


# In[ ]:


droprepeatmasker


# # RepeatMakser LINE TSNE and PCA

# In[ ]:


X_2u['class'] = np.where(droprepeatmasker['class'].str.contains('LINE'),'LINE', droprepeatmasker['class'])
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u, x="a", y="b", hue='class')
plt.title('Upstream TSNE with repearmasker class LINE')
plt.show()


# In[ ]:


X_2d['class'] = np.where(droprepeatmasker['class'].str.contains('LINE'),'LINE', droprepeatmasker['class'])
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='class')
plt.title('Downstream TSNE with repearmasker class LINE')
plt.show()


# # RepeatMakser Alu TSNE and PCA

# In[ ]:


# downstream
X_2d['class'] = np.where(droprepeatmasker['class'].str.contains('Alu'),'Alu', droprepeatmasker['class'])
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='class')
plt.title('Downstream TSNE with repearmasker class Alu')
plt.show()


# In[ ]:


# upstream
droprepeatmasker['superclass'] = np.where(droprepeatmasker['class'].str.contains('LINE'),'LINE', droprepeatmasker['superclass'])
droprepeatmasker['superclass'] = np.where(droprepeatmasker['class'].str.contains('Alu'),'Alu', droprepeatmasker['superclass'])
droprepeatmasker['superclass'] = np.where(droprepeatmasker['class'].str.contains('LTR'),'LTR', droprepeatmasker['superclass'])
droprepeatmasker['superclass'] = np.where(droprepeatmasker['class'].str.contains('MER'),'MER', droprepeatmasker['superclass'])
droprepeatmasker['superclass'] = np.where(droprepeatmasker['class'].str.contains('HERV'),'HERV', droprepeatmasker['superclass'])
X_2u['superclass'] = droprepeatmasker['superclass']

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u, x="a", y="b", hue='superclass')
plt.title('Upstream TSNE with repearmasker class Alu')
plt.show()


# In[ ]:


# downstream
droprepeatmasker['superclass'] = np.where(droprepeatmasker['class'].str.contains('LINE'),'LINE', droprepeatmasker['superclass'])
droprepeatmasker['superclass'] = np.where(droprepeatmasker['class'].str.contains('Alu'),'Alu', droprepeatmasker['superclass'])
droprepeatmasker['superclass'] = np.where(droprepeatmasker['class'].str.contains('LTR'),'LTR', droprepeatmasker['superclass'])
droprepeatmasker['superclass'] = np.where(droprepeatmasker['class'].str.contains('MER'),'MER', droprepeatmasker['superclass'])
droprepeatmasker['superclass'] = np.where(droprepeatmasker['class'].str.contains('HERV'),'HERV', droprepeatmasker['superclass'])
X_2d['superclass'] = droprepeatmasker['superclass']

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='superclass')
plt.title('Downstream TSNE with repearmasker class Alu')
plt.show()


# # RepeatMakser SubAlu Class PCA and TSNE

# In[ ]:


droprepeatmasker['subAluclass'] = np.where(droprepeatmasker['class'].str.contains('AluY'),'AluY', droprepeatmasker['subAluclass'])
droprepeatmasker['subAluclass'] = np.where(droprepeatmasker['class'].str.contains('AluS'),'AluS', droprepeatmasker['subAluclass'])
droprepeatmasker['subAluclass'] = np.where(droprepeatmasker['class'].str.contains('AluJ'),'AluJ', droprepeatmasker['subAluclass'])

X_2d['subAluclass'] = droprepeatmasker['subAluclass']

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='subAluclass')
plt.title('Downstream TSNE with repearmasker subclass Alu')
plt.show()


# In[ ]:


droprepeatmasker['subAluclass'] = np.where(droprepeatmasker['class'].str.contains('AluY'),'AluY', droprepeatmasker['subAluclass'])
droprepeatmasker['subAluclass'] = np.where(droprepeatmasker['class'].str.contains('AluS'),'AluS', droprepeatmasker['subAluclass'])
droprepeatmasker['subAluclass'] = np.where(droprepeatmasker['class'].str.contains('AluJ'),'AluJ', droprepeatmasker['subAluclass'])

X_2u['subAluclass'] = droprepeatmasker['subAluclass']

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u, x="a", y="b", hue='subAluclass')
plt.title('Upstream TSNE with repearmasker subclass Alu')
plt.show()


# In[ ]:


droprepeatmasker['subAluclass'] = np.where(droprepeatmasker['class'].str.contains('AluY'),'AluY', droprepeatmasker['subAluclass'])
droprepeatmasker['subAluclass'] = np.where(droprepeatmasker['class'].str.contains('AluS'),'AluS', droprepeatmasker['subAluclass'])
droprepeatmasker['subAluclass'] = np.where(droprepeatmasker['class'].str.contains('AluJ'),'AluJ', droprepeatmasker['subAluclass'])

principalupstreamwinnum['subAluclass'] = droprepeatmasker['subAluclass']

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principalupstreamwinnum, x="a", y="b", hue='subAluclass')
plt.title('Upstream PCA with repearmasker subclass Alu')
plt.show()


# In[ ]:


droprepeatmasker['subAluclass'] = np.where(droprepeatmasker['class'].str.contains('AluY'),'AluY', droprepeatmasker['subAluclass'])
droprepeatmasker['subAluclass'] = np.where(droprepeatmasker['class'].str.contains('AluS'),'AluS', droprepeatmasker['subAluclass'])
droprepeatmasker['subAluclass'] = np.where(droprepeatmasker['class'].str.contains('AluJ'),'AluJ', droprepeatmasker['subAluclass'])

principaldownstreamwinnum['subAluclass'] = droprepeatmasker['subAluclass']

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue='subAluclass')
plt.title('Downstream PCA with repearmasker subclass Alu')
plt.show()


# In[ ]:


droprepeatmasker['superclass']= None
droprepeatmasker


# In[ ]:


droprepeatmasker['subAluclass']= None
droprepeatmasker


# In[ ]:


# AluS, AluJ, AluY 
# flatten all the windows - instead of concat, sum/ mean  (cluster)
# 


# # chrXch38 Validation Coloring

# In[17]:


import numpy as np
import pandas as pd
ch38 = pd.read_csv("/data/projects/nanopore/RepeatExpansion/TR_downstreamAnalysis/chrXCh38/closestsortedchrXGRCh38.bed", sep="\t", header=None, names = ['chrX','1','2','chrX1','3','4','class','distance'])
ch38


# In[18]:


ch38drop=ch38.drop_duplicates(subset=["1","2"],keep="first") 
ch38drop


# In[39]:


# pca downstream 
c

# pca upstream
principaldownstreamwinnum['class'] =ch38drop['class']
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue='class')
plt.title("downstream pca with ch38")
plt.show()

#tsne down
X_2d['class']=np.where(ch38drop['class'].str.contains("P"),"P",ch38drop['class'])
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='class')
plt.title('Downstream TSNE with ch38')
plt.show()

#tsne up
X_2u['class']=np.where(ch38drop['class'].str.contains("P"),"P",ch38drop['class'])
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u, x="a", y="b", hue='class')
plt.title('Upstream TSNE with ch38')
plt.show()


# In[40]:


ch38drop['superclass']= None
ch38drop


# In[41]:


# upstream
ch38drop['superclass'] = np.where(ch38drop['class'].str.contains('pELS'),'pELS', ch38drop['superclass'])
ch38drop['superclass'] = np.where(ch38drop['class'].str.contains('CTCF-only'),'CTCF', ch38drop['superclass'])
ch38drop['superclass'] = np.where(ch38drop['class'].str.contains('P'),'PLS', ch38drop['superclass'])
ch38drop['superclass'] = np.where(ch38drop['class'].str.contains('DNase'),'DNase', ch38drop['superclass'])
X_2u['superclass'] = ch38drop['superclass']

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u, x="a", y="b", hue='superclass')
plt.title('Upstream TSNE with ch38drop class Alu')
plt.show()


# In[42]:


ch38drop['superclass'] = np.where(ch38drop['class'].str.contains('pELS'),'pELS', ch38drop['superclass'])
ch38drop['superclass'] = np.where(ch38drop['class'].str.contains('CTCF-only'),'CTCF', ch38drop['superclass'])
ch38drop['superclass'] = np.where(ch38drop['class'].str.contains('P'),'PLS', ch38drop['superclass'])
ch38drop['superclass'] = np.where(ch38drop['class'].str.contains('DNase'),'DNase', ch38drop['superclass'])
X_2d['superclass'] = ch38drop['superclass']

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='superclass')
plt.title('Downstream TSNE with ch38drop class Alu')
plt.show()


# In[43]:


ch38drop['superclass'] = np.where(ch38drop['class'].str.contains('CTCF-bound'),'CTCFbound', ch38drop['superclass'])
X_2u1['superclass'] = ch38drop['superclass']

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u1, x="a", y="b", hue='superclass')
plt.title('Downstream TSNE with ch38drop class CTCFbound')
plt.show()


# In[44]:


ch38drop


# In[ ]:


ch38drop.head(40)


# # ch38 CTCF BOUND Coloring 

# In[45]:


ch38drop['CTCFbound']= None
ch38drop


# In[19]:


#Upsteam tsne
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('DNase-H3K4me3,CTCF-bound'),'DNase-H3K4me3', ch38drop['CTCFbound'])
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('PLS,CTCF-bound'),'PLS', ch38drop['CTCFbound'])
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('pELS,CTCF-bound'),'pELS', ch38drop['CTCFbound'])
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('CTCF-only,CTCF-bound'),'CTCFonly', ch38drop['CTCFbound'])
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('dELS,CTCF-bound'),'dELS', ch38drop['CTCFbound'])
X_2u['CTCFbound'] = ch38drop['CTCFbound']

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u, x="a", y="b", hue='CTCFbound')
plt.title('Downstream TSNE with ch38drop class CTCFbound')
plt.show()


# In[20]:


#Downsteam tsne
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('DNase-H3K4me3,CTCF-bound'),'DNase-H3K4me3', ch38drop['CTCFbound'])
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('PLS,CTCF-bound'),'PLS', ch38drop['CTCFbound'])
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('pELS,CTCF-bound'),'pELS', ch38drop['CTCFbound'])
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('CTCF-only,CTCF-bound'),'CTCFonly', ch38drop['CTCFbound'])
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('dELS,CTCF-bound'),'dELS', ch38drop['CTCFbound'])
X_2d['CTCFbound'] = ch38drop['CTCFbound']

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='CTCFbound')
plt.title('Downstream TSNE with ch38drop class CTCFbound')
plt.show()


# In[46]:


#downsteam pca
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('DNase-H3K4me3,CTCF-bound'),'DNase-H3K4me3', ch38drop['CTCFbound'])
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('PLS,CTCF-bound'),'PLS', ch38drop['CTCFbound'])
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('pELS,CTCF-bound'),'pELS', ch38drop['CTCFbound'])
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('CTCF-only,CTCF-bound'),'CTCFonly', ch38drop['CTCFbound'])
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('dELS,CTCF-bound'),'dELS', ch38drop['CTCFbound'])
principaldownstreamwinnum['CTCFbound'] = ch38drop['CTCFbound']

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue='CTCFbound')
plt.title('Downstream principaldownstreamwinnum with ch38drop class CTCFbound')
plt.show()


# In[22]:


#upsteam pca
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('DNase-H3K4me3,CTCF-bound'),'DNase-H3K4me3', ch38drop['CTCFbound'])
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('PLS,CTCF-bound'),'PLS', ch38drop['CTCFbound'])
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('pELS,CTCF-bound'),'pELS', ch38drop['CTCFbound'])
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('CTCF-only,CTCF-bound'),'CTCFonly', ch38drop['CTCFbound'])
ch38drop['CTCFbound'] = np.where(ch38drop['class'].str.contains('dELS,CTCF-bound'),'dELS', ch38drop['CTCFbound'])
principalupstreamwinnum['CTCFbound'] = ch38drop['CTCFbound']

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principalupstreamwinnum, x="a", y="b", hue='CTCFbound')
plt.title('Downstream principalupstreamwinnum with ch38drop class CTCFbound')
plt.show()


# # filter by the distance (ch38)

# In[39]:


X_2u['distance'] = np.where((ch38drop['distance']>-2000)&(ch38drop['distance']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u, x="a", y="b", hue='distance')
plt.title('Upstream TSNE with ch38drop distance of larger than -2000')
plt.show()


# In[34]:


X_2u1['distance'] = np.where((ch38drop['distance']>-2000)&(ch38drop['distance']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u1, x="a", y="b", hue='distance')
plt.title('Upstream TSNE with ch38drop distance of larger than -2000')
plt.show()


# In[50]:


X_2u['distance'] = np.where((ch38drop['distance']>-2000)&(ch38drop['distance']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2u, x="a", y="b", hue='distance')
plt.title('Upstream TSNE with ch38drop distance of larger than -2000')
plt.show()


# In[37]:


principalupstreamwinnum['distance'] = np.where((ch38drop['distance']>-2000)&(ch38drop['distance']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principalupstreamwinnum, x="a", y="b", hue='distance')
plt.title('Downstream with ch38drop distance of larger than -2000')
plt.show()


# In[53]:


principalupstreamwinnum['distance'] = np.where((ch38drop['distance']>-1000)&(ch38drop['distance']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principalupstreamwinnum, x="a", y="b", hue='distance')
plt.title('Downstream with ch38drop distance of larger than -2000')
plt.show()


# In[55]:


principalupstreamwinnum['distance'] = np.where((ch38drop['distance']>-500)&(ch38drop['distance']<-0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principalupstreamwinnum, x="a", y="b", hue='distance')
plt.title('Downstream with ch38drop distance of larger than -2000')
plt.show()


# # validation Downstream (filter by the distance)

# In[47]:


# 2000
principaldownstreamwinnum['distance'] = np.where((ch38drop['distance']<2000)&(ch38drop['distance']>0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue='distance')
plt.title('Downstream with ch38drop distance of larger than -2000')
plt.show()


# In[52]:


# 50000
principaldownstreamwinnum['distance'] = np.where((ch38drop['distance']<50000)&(ch38drop['distance']>0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue='distance')
plt.title('Downstream with ch38drop distance of smaller than 5000')
plt.show()


# In[48]:


# 5000
principaldownstreamwinnum['distance'] = np.where((ch38drop['distance']<5000)&(ch38drop['distance']>0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue='distance')
plt.title('Downstream with ch38drop distance of smaller than 5000')
plt.show()


# In[49]:


# 1000
principaldownstreamwinnum['distance'] = np.where((ch38drop['distance']<1000)&(ch38drop['distance']>0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue='distance')
plt.title('Downstream with ch38drop distance of smaller than 1000')
plt.show()


# In[50]:


# 500
principaldownstreamwinnum['distance'] = np.where((ch38drop['distance']<500)&(ch38drop['distance']>0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnum, x="a", y="b", hue='distance')
plt.title('Downstream with ch38drop distance of smaller than 500')
plt.show()


# In[51]:


X_2d['distance'] = np.where((ch38drop['distance']<2000)&(ch38drop['distance']>0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='distance')
plt.title('Downstream TSNE with ch38drop larger than -2000')
plt.show()


# In[51]:


X_2d['distance'] = np.where((ch38drop['distance']<1000)&(ch38drop['distance']>0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='distance')
plt.title('Downstream TSNE with ch38drop smaller than 1000')
plt.show()pc


# In[58]:


X_2d['distance'] = np.where((ch38drop['distance']<500)&(ch38drop['distance']>0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=X_2d, x="a", y="b", hue='distance')
plt.title('Downstream TSNE with ch38drop smaller than 500')
plt.show()


# # UMAP Downstream

# In[72]:


import numpy as np
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
get_ipython().run_line_magic('matplotlib', 'inline')


# In[25]:


get_ipython().system('conda install -c conda-forge umap-learn -y')


# In[73]:


# Data manipulation
import pandas as pd # for data manipulation
import numpy as np # for data manipulation

# Visualization
import plotly.express as px # for data visualization
import matplotlib.pyplot as plt # for showing handwritten digits

# Skleran
from sklearn.datasets import load_digits # for MNIST data
from sklearn.model_selection import train_test_split # for splitting data into train and test samples

# UMAP dimensionality reduction
from umap import UMAP

import umap
# In[27]:


get_ipython().system('pip install umap-learn')


# In[14]:


Downstreamdropwin


# In[15]:


Upstreamdropwin


# In[21]:


Downstreamdropwin
x = Downstreamdropwin.values
y = Downstreamdropwin.values
x = StandardScaler().fit_transform(x)

downstreamumap = UMAP(n_neighbors=100, # default 15, The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation.
               n_components=2, # default 2, The dimension of the space to embed into.
               metric='euclidean', # default 'euclidean', The metric to use to compute distances in high dimensional space.
               n_epochs=1000, # default None, The number of training epochs to be used in optimizing the low dimensional embedding. Larger values result in more accurate embeddings. 
               learning_rate=1.0, # default 1.0, The initial learning rate for the embedding optimization.
               init='spectral', # default 'spectral', How to initialize the low dimensional embedding. Options are: {'spectral', 'random', A numpy array of initial embedding positions}.
               min_dist=0.1, # default 0.1, The effective minimum distance between embedded points.
               spread=1.0, # default 1.0, The effective scale of embedded points. In combination with ``min_dist`` this determines how clustered/clumped the embedded points are.
               low_memory=False, # default False, For some datasets the nearest neighbor computation can consume a lot of memory. If you find that UMAP is failing due to memory constraints consider setting this option to True.
               set_op_mix_ratio=1.0, # default 1.0, The value of this parameter should be between 0.0 and 1.0; a value of 1.0 will use a pure fuzzy union, while 0.0 will use a pure fuzzy intersection.
               local_connectivity=1, # default 1, The local connectivity required -- i.e. the number of nearest neighbors that should be assumed to be connected at a local level.
               repulsion_strength=1.0, # default 1.0, Weighting applied to negative samples in low dimensional embedding optimization.
               negative_sample_rate=5, # default 5, Increasing this value will result in greater repulsive force being applied, greater optimization cost, but slightly more accuracy.
               transform_queue_size=4.0, # default 4.0, Larger values will result in slower performance but more accurate nearest neighbor evaluation.
               a=None, # default None, More specific parameters controlling the embedding. If None these values are set automatically as determined by ``min_dist`` and ``spread``.
               b=None, # default None, More specific parameters controlling the embedding. If None these values are set automatically as determined by ``min_dist`` and ``spread``.
               random_state=42, # default: None, If int, random_state is the seed used by the random number generator;
               metric_kwds=None, # default None) Arguments to pass on to the metric, such as the ``p`` value for Minkowski distance.
               angular_rp_forest=False, # default False, Whether to use an angular random projection forest to initialise the approximate nearest neighbor search.
               target_n_neighbors=-1, # default -1, The number of nearest neighbors to use to construct the target simplcial set. If set to -1 use the ``n_neighbors`` value.
               #target_metric='categorical', # default 'categorical', The metric used to measure distance for a target array is using supervised dimension reduction. By default this is 'categorical' which will measure distance in terms of whether categories match or are different. 
               #target_metric_kwds=None, # dict, default None, Keyword argument to pass to the target metric when performing supervised dimension reduction. If None then no arguments are passed on.
               #target_weight=0.5, # default 0.5, weighting factor between data topology and target topology.
               transform_seed=42, # default 42, Random seed used for the stochastic aspects of the transform operation.
               verbose=False, # default False, Controls verbosity of logging.
               unique=False, # default False, Controls if the rows of your data should be uniqued before being embedded. 
              )

# Fit and transform the data
Xdown = downstreamumap.fit_transform(x)

# Check the shape of the new data
print('Shape of X_trans: ', Xdown.shape)



# In[17]:


Xdowndf=pd.DataFrame(Xdown)


# In[18]:


plt.scatter(x=Xdowndf[0],y=Xdowndf[1])
plt.show()


# In[74]:



## hue corresponds to the value
## write the function that would do everything 
## dark - most similar features 
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=Xdowndf,x=Xdowndf[0],y=Xdowndf[1], hue=Xdowndf[0])


# In[ ]:


for key in translater.keys():

    Xdowndf[key] = sumTR[key]
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.scatterplot(data=Xdowndf, x=Xdowndf[0],y=Xdowndf[1], hue=key)
    plt.title(key)
    plt.show()
    


# In[19]:


Xdowndf.astype(float)


# In[20]:


plt.hist(Xdowndf.astype(float))


# In[41]:


from sklearn.datasets import fetch_openml
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

# Dimension reduction and clustering libraries
import umap
import sklearn.cluster as cluster
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score


# # UMAP upstream

# In[74]:


Downstreamdropwin
x = Upstreamdropwin.values
y = Upstreamdropwin.values
x = StandardScaler().fit_transform(x)

Upstreamumap = UMAP(n_neighbors=100, # default 15, The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation.
               n_components=2, # default 2, The dimension of the space to embed into.
               metric='euclidean', # default 'euclidean', The metric to use to compute distances in high dimensional space.
               n_epochs=1000, # default None, The number of training epochs to be used in optimizing the low dimensional embedding. Larger values result in more accurate embeddings. 
               learning_rate=1.0, # default 1.0, The initial learning rate for the embedding optimization.
               init='spectral', # default 'spectral', How to initialize the low dimensional embedding. Options are: {'spectral', 'random', A numpy array of initial embedding positions}.
               min_dist=0.1, # default 0.1, The effective minimum distance between embedded points.
               spread=1.0, # default 1.0, The effective scale of embedded points. In combination with ``min_dist`` this determines how clustered/clumped the embedded points are.
               low_memory=False, # default False, For some datasets the nearest neighbor computation can consume a lot of memory. If you find that UMAP is failing due to memory constraints consider setting this option to True.
               set_op_mix_ratio=1.0, # default 1.0, The value of this parameter should be between 0.0 and 1.0; a value of 1.0 will use a pure fuzzy union, while 0.0 will use a pure fuzzy intersection.
               local_connectivity=1, # default 1, The local connectivity required -- i.e. the number of nearest neighbors that should be assumed to be connected at a local level.
               repulsion_strength=1.0, # default 1.0, Weighting applied to negative samples in low dimensional embedding optimization.
               negative_sample_rate=5, # default 5, Increasing this value will result in greater repulsive force being applied, greater optimization cost, but slightly more accuracy.
               transform_queue_size=4.0, # default 4.0, Larger values will result in slower performance but more accurate nearest neighbor evaluation.
               a=None, # default None, More specific parameters controlling the embedding. If None these values are set automatically as determined by ``min_dist`` and ``spread``.
               b=None, # default None, More specific parameters controlling the embedding. If None these values are set automatically as determined by ``min_dist`` and ``spread``.
               random_state=42, # default: None, If int, random_state is the seed used by the random number generator;
               metric_kwds=None, # default None) Arguments to pass on to the metric, such as the ``p`` value for Minkowski distance.
               angular_rp_forest=False, # default False, Whether to use an angular random projection forest to initialise the approximate nearest neighbor search.
               target_n_neighbors=-1, # default -1, The number of nearest neighbors to use to construct the target simplcial set. If set to -1 use the ``n_neighbors`` value.
               #target_metric='categorical', # default 'categorical', The metric used to measure distance for a target array is using supervised dimension reduction. By default this is 'categorical' which will measure distance in terms of whether categories match or are different. 
               #target_metric_kwds=None, # dict, default None, Keyword argument to pass to the target metric when performing supervised dimension reduction. If None then no arguments are passed on.
               #target_weight=0.5, # default 0.5, weighting factor between data topology and target topology.
               transform_seed=42, # default 42, Random seed used for the stochastic aspects of the transform operation.
               verbose=False, # default False, Controls verbosity of logging.
               unique=False, # default False, Controls if the rows of your data should be uniqued before being embedded. 
              )

# Fit and transform the data
XUp = Upstreamumap.fit_transform(x)

# Check the shape of the new data
print('Shape of X_trans: ', XUp.shape)


# In[75]:


XUpdf=pd.DataFrame(XUp)


# In[76]:


for key in translater.keys():

    XUpdf[key] = sumTRUP[key]
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.scatterplot(data=XUpdf, x=XUpdf[0],y=XUpdf[1], hue=key)
    plt.title(key)
    plt.show()


# # Louvain clustering

# In[ ]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created in November 2018
@author: Nathan de Lara <nathan.delara@polytechnique.org>
@author: Quentin Lutz <qlutz@enst.fr>
@author: Thomas Bonald <bonald@enst.fr>
"""
from typing import Union, Optional

import numpy as np
from scipy import sparse

from sknetwork.clustering.base import BaseClustering
from sknetwork.clustering.louvain_core import fit_core
from sknetwork.clustering.postprocess import reindex_labels
from sknetwork.utils.check import check_random_state, get_probs
from sknetwork.utils.format import check_format, get_adjacency, directed2undirected
from sknetwork.utils.membership import get_membership
from sknetwork.utils.verbose import VerboseMixin


[docs]class Louvain(BaseClustering, VerboseMixin):
    """Louvain algorithm for clustering graphs by maximization of modularity.

    For bipartite graphs, the algorithm maximizes Barber's modularity by default.

    Parameters
    ----------
    resolution :
        Resolution parameter.
    modularity : str
        Which objective function to maximize. Can be ``'Dugue'``, ``'Newman'`` or ``'Potts'`` (default = ``'dugue'``).
    tol_optimization :
        Minimum increase in the objective function to enter a new optimization pass.
    tol_aggregation :
        Minimum increase in the objective function to enter a new aggregation pass.
    n_aggregations :
        Maximum number of aggregations.
        A negative value is interpreted as no limit.
    shuffle_nodes :
        Enables node shuffling before optimization.
    sort_clusters :
        If ``True``, sort labels in decreasing order of cluster size.
    return_probs :
        If ``True``, return the probability distribution over clusters (soft clustering).
    return_aggregate :
        If ``True``, return the adjacency matrix of the graph between clusters.
    random_state :
        Random number generator or random seed. If None, numpy.random is used.
    verbose :
        Verbose mode.

    Attributes
    ----------
    labels_ : np.ndarray, shape (n_labels,)
        Label of each node.
    probs_ : sparse.csr_matrix, shape (n_row, n_labels)
        Probability distribution over labels.
    labels_row_, labels_col_ : np.ndarray
        Labels of rows and columns, for bipartite graphs.
    probs_row_, probs_col_ : sparse.csr_matrix, shape (n_row, n_labels)
        Probability distributions over labels for rows and columns (for bipartite graphs).
    aggregate_ : sparse.csr_matrix
        Aggregate adjacency matrix or biadjacency matrix between clusters.

    Example
    -------
    >>> from sknetwork.clustering import Louvain
    >>> from sknetwork.data import karate_club
    >>> louvain = Louvain()
    >>> adjacency = karate_club()
    >>> labels = louvain.fit_predict(adjacency)
    >>> len(set(labels))
    4

    References
    ----------
    * Blondel, V. D., Guillaume, J. L., Lambiotte, R., & Lefebvre, E. (2008).
      `Fast unfolding of communities in large networks.
      <https://arxiv.org/abs/0803.0476>`_
      Journal of statistical mechanics: theory and experiment, 2008.

    * Dugué, N., & Perez, A. (2015).
      `Directed Louvain: maximizing modularity in directed networks
      <https://hal.archives-ouvertes.fr/hal-01231784/document>`_
      (Doctoral dissertation, Université d'Orléans).

    * Barber, M. J. (2007).
      `Modularity and community detection in bipartite networks
      <https://arxiv.org/pdf/0707.1616>`_
      Physical Review E, 76(6).
    """
    def __init__(self, resolution: float = 1, modularity: str = 'dugue', tol_optimization: float = 1e-3,
                 tol_aggregation: float = 1e-3, n_aggregations: int = -1, shuffle_nodes: bool = False,
                 sort_clusters: bool = True, return_probs: bool = True, return_aggregate: bool = True,
                 random_state: Optional[Union[np.random.RandomState, int]] = None, verbose: bool = False):
        super(Louvain, self).__init__(sort_clusters=sort_clusters, return_probs=return_probs,
                                      return_aggregate=return_aggregate)
        VerboseMixin.__init__(self, verbose)

        self.labels_ = None
        self.resolution = resolution
        self.modularity = modularity.lower()
        self.tol = tol_optimization
        self.tol_aggregation = tol_aggregation
        self.n_aggregations = n_aggregations
        self.shuffle_nodes = shuffle_nodes
        self.random_state = check_random_state(random_state)
        self.bipartite = None

    def _optimize(self, adjacency_norm, probs_ou, probs_in):
        """One local optimization pass of the Louvain algorithm

        Parameters
        ----------
        adjacency_norm :
            the norm of the adjacency
        probs_ou :
            the array of degrees of the adjacency
        probs_in :
            the array of degrees of the transpose of the adjacency

        Returns
        -------
        labels :
            the communities of each node after optimization
        pass_increase :
            the increase in modularity gained after optimization
        """
        node_probs_in = probs_in.astype(np.float32)
        node_probs_ou = probs_ou.astype(np.float32)

        adjacency = 0.5 * directed2undirected(adjacency_norm)

        self_loops = adjacency.diagonal().astype(np.float32)

        indptr: np.ndarray = adjacency.indptr
        indices: np.ndarray = adjacency.indices
        data: np.ndarray = adjacency.data.astype(np.float32)

        return fit_core(self.resolution, self.tol, node_probs_ou, node_probs_in, self_loops, data, indices, indptr)

    @staticmethod
    def _aggregate(adjacency_norm, probs_out, probs_in, membership: Union[sparse.csr_matrix, np.ndarray]):
        """Aggregate nodes belonging to the same cluster.

        Parameters
        ----------
        adjacency_norm :
            the norm of the adjacency
        probs_out :
            the array of degrees of the adjacency
        probs_in :
            the array of degrees of the transpose of the adjacency
        membership :
            membership matrix (rows).

        Returns
        -------
        Aggregate graph.
        """
        adjacency_norm = (membership.T.dot(adjacency_norm.dot(membership))).tocsr()
        probs_in = np.array(membership.T.dot(probs_in).T)
        probs_out = np.array(membership.T.dot(probs_out).T)
        return adjacency_norm, probs_out, probs_in

[docs]    def fit(self, input_matrix: Union[sparse.csr_matrix, np.ndarray], force_bipartite: bool = False) -> 'Louvain':
        """Fit algorithm to data.

        Parameters
        ----------
        input_matrix :
            Adjacency matrix or biadjacency matrix of the graph.
        force_bipartite :
            If ``True``, force the input matrix to be considered as a biadjacency matrix even if square.

        Returns
        -------
        self: :class:`Louvain`
        """
        self._init_vars()
        input_matrix = check_format(input_matrix)
        if self.modularity == 'dugue':
            adjacency, self.bipartite = get_adjacency(input_matrix, force_directed=True,
                                                      force_bipartite=force_bipartite)
        else:
            adjacency, self.bipartite = get_adjacency(input_matrix, force_bipartite=force_bipartite)

        n = adjacency.shape[0]

        index = np.arange(n)
        if self.shuffle_nodes:
            index = self.random_state.permutation(index)
            adjacency = adjacency[index][:, index]

        if self.modularity == 'potts':
            probs_out = get_probs('uniform', adjacency)
            probs_in = probs_out.copy()
        elif self.modularity == 'newman':
            probs_out = get_probs('degree', adjacency)
            probs_in = probs_out.copy()
        elif self.modularity == 'dugue':
            probs_out = get_probs('degree', adjacency)
            probs_in = get_probs('degree', adjacency.T)
        else:
            raise ValueError('Unknown modularity function.')

        adjacency_cluster = adjacency / adjacency.data.sum()

        membership = sparse.identity(n, format='csr')
        increase = True
        count_aggregations = 0
        self.log.print("Starting with", n, "nodes.")
        while increase:
            count_aggregations += 1

            labels_cluster, pass_increase = self._optimize(adjacency_cluster, probs_out, probs_in)
            _, labels_cluster = np.unique(labels_cluster, return_inverse=True)

            if pass_increase <= self.tol_aggregation:
                increase = False
            else:
                membership_cluster = get_membership(labels_cluster)
                membership = membership.dot(membership_cluster)
                adjacency_cluster, probs_out, probs_in = self._aggregate(adjacency_cluster, probs_out, probs_in,
                                                                         membership_cluster)

                n = adjacency_cluster.shape[0]
                if n == 1:
                    break
            self.log.print("Aggregation", count_aggregations, "completed with", n, "clusters and ",
                           pass_increase, "increment.")
            if count_aggregations == self.n_aggregations:
                break

        if self.sort_clusters:
            labels = reindex_labels(membership.indices)
        else:
            labels = membership.indices
        if self.shuffle_nodes:
            reverse = np.empty(index.size, index.dtype)
            reverse[index] = np.arange(index.size)
            labels = labels[reverse]

        self.labels_ = labels
        if self.bipartite:
            self._split_vars(input_matrix.shape)
        self._secondary_outputs(input_matrix)

        return self


# # Plt Sum of SEQ TR DownStream

# In[106]:


dfDownwinnum


# In[107]:


def sum_df(df, tr_id, seqclass):
    #print (df)
    dfDownwinnum = df.loc[df["Result"] == (tr_id)].copy()
    sumTRpltdown =dfDownwinnum[str(seqclass)].sum()
    return sumTRpltdown



colsDown = list()
colsDown = dfDownwinnum.columns.tolist()
print(type(colsDown))
print(len(colsDown))
#cols_newDOWN.extend(colsDOWN[0:]) 

#print(DownstreamMatrix.columns)
sumTRpltdown= {}

cnt = 0

for tr_id in list(set(dfDownwinnum["Result"].tolist())):
    sumTRpltdown[tr_id]=[]
    for seqclass in range(9,69):
        
        sumTRpltdown[tr_id].append(sum_df(dfDownwinnum,tr_id, seqclass))
    dfDownwinnumdrop = dfDownwinnum


# In[108]:


for tr_id in range(9,69):
    plt.plot(sumTRpltdown[tr_id])
    plt.title(tr_id)
    plt.show()


# In[109]:


dfsumTRpltdown=pd.DataFrame(sumTRpltdown).T


# In[110]:


dfsumTRpltdown


# In[35]:



x = dfsumTRpltdown.values
y = dfsumTRpltdown.values
x = StandardScaler().fit_transform(x)

dfsumTRpltdownumap = UMAP(n_neighbors=100, # default 15, The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation.
               n_components=2, # default 2, The dimension of the space to embed into.
               metric='euclidean', # default 'euclidean', The metric to use to compute distances in high dimensional space.
               n_epochs=1000, # default None, The number of training epochs to be used in optimizing the low dimensional embedding. Larger values result in more accurate embeddings. 
               learning_rate=1.0, # default 1.0, The initial learning rate for the embedding optimization.
               init='spectral', # default 'spectral', How to initialize the low dimensional embedding. Options are: {'spectral', 'random', A numpy array of initial embedding positions}.
               min_dist=0.1, # default 0.1, The effective minimum distance between embedded points.
               spread=1.0, # default 1.0, The effective scale of embedded points. In combination with ``min_dist`` this determines how clustered/clumped the embedded points are.
               low_memory=False, # default False, For some datasets the nearest neighbor computation can consume a lot of memory. If you find that UMAP is failing due to memory constraints consider setting this option to True.
               set_op_mix_ratio=1.0, # default 1.0, The value of this parameter should be between 0.0 and 1.0; a value of 1.0 will use a pure fuzzy union, while 0.0 will use a pure fuzzy intersection.
               local_connectivity=1, # default 1, The local connectivity required -- i.e. the number of nearest neighbors that should be assumed to be connected at a local level.
               repulsion_strength=1.0, # default 1.0, Weighting applied to negative samples in low dimensional embedding optimization.
               negative_sample_rate=5, # default 5, Increasing this value will result in greater repulsive force being applied, greater optimization cost, but slightly more accuracy.
               transform_queue_size=4.0, # default 4.0, Larger values will result in slower performance but more accurate nearest neighbor evaluation.
               a=None, # default None, More specific parameters controlling the embedding. If None these values are set automatically as determined by ``min_dist`` and ``spread``.
               b=None, # default None, More specific parameters controlling the embedding. If None these values are set automatically as determined by ``min_dist`` and ``spread``.
               random_state=42, # default: None, If int, random_state is the seed used by the random number generator;
               metric_kwds=None, # default None) Arguments to pass on to the metric, such as the ``p`` value for Minkowski distance.
               angular_rp_forest=False, # default False, Whether to use an angular random projection forest to initialise the approximate nearest neighbor search.
               target_n_neighbors=-1, # default -1, The number of nearest neighbors to use to construct the target simplcial set. If set to -1 use the ``n_neighbors`` value.
               #target_metric='categorical', # default 'categorical', The metric used to measure distance for a target array is using supervised dimension reduction. By default this is 'categorical' which will measure distance in terms of whether categories match or are different. 
               #target_metric_kwds=None, # dict, default None, Keyword argument to pass to the target metric when performing supervised dimension reduction. If None then no arguments are passed on.
               #target_weight=0.5, # default 0.5, weighting factor between data topology and target topology.
               transform_seed=42, # default 42, Random seed used for the stochastic aspects of the transform operation.
               verbose=False, # default False, Controls verbosity of logging.
               unique=False, # default False, Controls if the rows of your data should be uniqued before being embedded. 
              )

# Fit and transform the data
Xdownumap = dfsumTRpltdownumap.fit_transform(x)

# Check the shape of the new data
print('Shape of X_trans: ', Xdownumap.shape)


# In[36]:


Xdownumapdf=pd.DataFrame(Xdownumap)


# In[77]:


# 
for key in translater.keys():

    Xdownumapdf[key] = sumTR[key]
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.scatterplot(data=Xdownumapdf, x=Xdownumapdf[0],y=Xdownumapdf[1], hue=key)
    plt.title(key)
    plt.show()


# # Plt Sum of SEQ TR UpStream

# In[98]:


def sum_df(df, tr_id, seqclass):
    #print (df)
    dfUpwinnum = df.loc[df["Result"] == (tr_id)].copy()
    sumTRpltup =dfUpwinnum[str(seqclass)].sum()
    return sumTRpltup


colsUp = list()
colsUp = dfUpwinnum.columns.tolist()
print(type(colsUp))
print(len(colsUp))
#cols_newDOWN.extend(colsDOWN[0:]) 

#print(DownstreamMatrix.columns)
sumTRpltup = {}

cnt = 0

for tr_id in list(set(dfUpwinnum["Result"].tolist())):
    sumTRpltup[tr_id]=[]
    for seqclass in range(9,69):
        
        sumTRpltup[tr_id].append(sum_df(dfUpwinnum,tr_id, seqclass))
    dfUpwinnumdrop = dfUpwinnum
    


# In[99]:


sumTRpltup


# In[ ]:



x = dfsumTRpltup.values
y = dfsumTRpltup.values
x = StandardScaler().fit_transform(x)

dfsumTRpltupumap = UMAP(n_neighbors=100, # default 15, The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation.
               n_components=2, # default 2, The dimension of the space to embed into.
               metric='euclidean', # default 'euclidean', The metric to use to compute distances in high dimensional space.
               n_epochs=1000, # default None, The number of training epochs to be used in optimizing the low dimensional embedding. Larger values result in more accurate embeddings. 
               learning_rate=1.0, # default 1.0, The initial learning rate for the embedding optimization.
               init='spectral', # default 'spectral', How to initialize the low dimensional embedding. Options are: {'spectral', 'random', A numpy array of initial embedding positions}.
               min_dist=0.1, # default 0.1, The effective minimum distance between embedded points.
               spread=1.0, # default 1.0, The effective scale of embedded points. In combination with ``min_dist`` this determines how clustered/clumped the embedded points are.
               low_memory=False, # default False, For some datasets the nearest neighbor computation can consume a lot of memory. If you find that UMAP is failing due to memory constraints consider setting this option to True.
               set_op_mix_ratio=1.0, # default 1.0, The value of this parameter should be between 0.0 and 1.0; a value of 1.0 will use a pure fuzzy union, while 0.0 will use a pure fuzzy intersection.
               local_connectivity=1, # default 1, The local connectivity required -- i.e. the number of nearest neighbors that should be assumed to be connected at a local level.
               repulsion_strength=1.0, # default 1.0, Weighting applied to negative samples in low dimensional embedding optimization.
               negative_sample_rate=5, # default 5, Increasing this value will result in greater repulsive force being applied, greater optimization cost, but slightly more accuracy.
               transform_queue_size=4.0, # default 4.0, Larger values will result in slower performance but more accurate nearest neighbor evaluation.
               a=None, # default None, More specific parameters controlling the embedding. If None these values are set automatically as determined by ``min_dist`` and ``spread``.
               b=None, # default None, More specific parameters controlling the embedding. If None these values are set automatically as determined by ``min_dist`` and ``spread``.
               random_state=42, # default: None, If int, random_state is the seed used by the random number generator;
               metric_kwds=None, # default None) Arguments to pass on to the metric, such as the ``p`` value for Minkowski distance.
               angular_rp_forest=False, # default False, Whether to use an angular random projection forest to initialise the approximate nearest neighbor search.
               target_n_neighbors=-1, # default -1, The number of nearest neighbors to use to construct the target simplcial set. If set to -1 use the ``n_neighbors`` value.
               #target_metric='categorical', # default 'categorical', The metric used to measure distance for a target array is using supervised dimension reduction. By default this is 'categorical' which will measure distance in terms of whether categories match or are different. 
               #target_metric_kwds=None, # dict, default None, Keyword argument to pass to the target metric when performing supervised dimension reduction. If None then no arguments are passed on.
               #target_weight=0.5, # default 0.5, weighting factor between data topology and target topology.
               transform_seed=42, # default 42, Random seed used for the stochastic aspects of the transform operation.
               verbose=False, # default False, Controls verbosity of logging.
               unique=False, # default False, Controls if the rows of your data should be uniqued before being embedded. 
              )

# Fit and transform the data
Xupumap = dfsumTRpltupumap.fit_transform(x)

# Check the shape of the new data
print('Shape of X_trans: ', Xupumap.shape)


# In[ ]:


dfsumTRpltup=pd.DataFrame(sumTRpltup).T


# In[57]:


Xupumapdf=pd.DataFrame(Xupumap)


# In[58]:


# 
for key in translater.keys():

    Xupumapdf[key] = sumTRUP[key]
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.scatterplot(data=Xupumapdf, x=Xupumapdf[0],y=Xupumapdf[1], hue=key)
    plt.title(key)
    plt.show()
    


# # Validation Data UMAP (Downstream)

# In[77]:


Xdownumapdf['distance'] = np.where((ch38drop['distance']<500)&(ch38drop['distance']>0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=Xdownumapdf, x=Xdownumapdf[0],y=Xdownumapdf[1], hue='distance')
plt.title('Downstream UMAP with ch38drop smaller than 500')
plt.show()


# In[78]:


Xdownumapdf['distance'] = np.where((ch38drop['distance']<5000)&(ch38drop['distance']>0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=Xdownumapdf, x=Xdownumapdf[0],y=Xdownumapdf[1], hue='distance')
plt.title('Downstream UMAP with ch38drop smaller than 5000')
plt.show()


# In[79]:


Xdownumapdf['distance'] = np.where((ch38drop['distance']<100)&(ch38drop['distance']>0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=Xdownumapdf, x=Xdownumapdf[0],y=Xdownumapdf[1], hue='distance')
plt.title('Downstream UMAP with ch38drop smaller than 100')
plt.show()


# In[80]:


Xdownumapdf['distance'] = np.where((ch38drop['distance']<1000)&(ch38drop['distance']>0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=Xdownumapdf, x=Xdownumapdf[0],y=Xdownumapdf[1], hue='distance')
plt.title('Downstream UMAP with ch38drop smaller than 1000')
plt.show()


# In[82]:


Xdownumapdf['distance'] = np.where((ch38drop['distance']<2000)&(ch38drop['distance']>0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=Xdownumapdf, x=Xdownumapdf[0],y=Xdownumapdf[1], hue='distance')
plt.title('Downstream UMAP with ch38drop smaller than 2000')
plt.show()


# # Validation Data UMAP (Upstream)

# In[73]:


Xupumapdf['distance'] = np.where((ch38drop['distance']>-500)&(ch38drop['distance']<0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=Xupumapdf, x=Xupumapdf[0],y=Xupumapdf[1], hue='distance')
plt.title('Upstream UMAP with ch38drop larger than -500')
plt.show()


# In[75]:


Xupumapdf['distance'] = np.where((ch38drop['distance']>-100)&(ch38drop['distance']<0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=Xupumapdf, x=Xupumapdf[0],y=Xupumapdf[1], hue='distance')
plt.title('Upstream UMAP with ch38drop larger than -100')
plt.show()


# In[74]:


Xupumapdf['distance'] = np.where((ch38drop['distance']>-5000)&(ch38drop['distance']<0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=Xupumapdf, x=Xupumapdf[0],y=Xupumapdf[1], hue='distance')
plt.title('Upstream UMAP with ch38drop larger than -5000')
plt.show()


# In[76]:


Xupumapdf['distance'] = np.where((ch38drop['distance']>-2000)&(ch38drop['distance']<0), True, False)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=Xupumapdf, x=Xupumapdf[0],y=Xupumapdf[1], hue='distance')
plt.title('Upstream UMAP with ch38drop larger than -2000')
plt.show()


#  tsne umap and pca are dimensionality reduction techniques, so they are used to project high dimensional data to lower dimensions
# which is helpful for visualization
# 
# clustering algorithms computationally group data points by their similarity, so the dimensionality reduced data would be an input to a clustering algorithm
# the reason we want to use a clustering algorithm (or community detection, they are kinda used interchangeably) is to identfy groups in our data that we may be able to see visually

# # Kmean 

# In[44]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


# # Kmean Upstream Matirx

# In[13]:


Upstreamdropwinscale=StandardScaler().fit_transform(Upstreamdropwin)


# In[14]:


Upstreamdropwinscale=StandardScaler().fit_transform(Upstreamdropwin)
#initialize kmeans parameters
kmeans_kwargs = {
"init": "random",
"n_init": 10,
"random_state": 1,
}

#create list to hold SSE values for each k
sse = []
for k in range(1, 11):
    kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
    kmeans.fit(Upstreamdropwinscale)
    sse.append(kmeans.inertia_)

#visualize results
plt.plot(range(1, 11), sse)
plt.xticks(range(1, 11))
plt.xlabel("Number of Clusters")
plt.ylabel("SSE")
plt.show()


# In[47]:


#instantiate the k-means class, using optimal number of clusters
kmeans = KMeans(init="random", n_clusters=3, n_init=10, random_state=1)

#fit k-means algorithm to data
kmeans.fit(Upstreamdropwinscale)

#view cluster assignments for each observation
kmeans.labels_


# In[50]:


upscaledff=pd.DataFrame(Upstreamdropwinscale)


# In[51]:


upscaledff


# In[52]:


#append cluster assingments to original DataFrame
upscaledff['cluster'] = kmeans.labels_

#view updated DataFrame
print(upscaledff)


# In[36]:


upscaledff['cluster'].unique()


# # Kmean Matrix coloring 

# In[ ]:


for key in translater.keys():

    upscaledff[key] = sumTRUP[key]
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.scatterplot(data=upscaledff, x="a", y="b", hue=key)
    plt.title(key)
    plt.show()


# In[42]:


sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=upscaledff, x=upscaledff['0'], y=upscaledff['1'], hue='cluster')
plt.title("matrix, kmean")
plt.show()


# # Kmean Upstream PCA

# In[55]:


# on PCA
principalupstreamwinnumscale=StandardScaler().fit_transform(principalupstreamwinnum)
#initialize kmeans parameters
kmeans_kwargs = {
"init": "random",
"n_init": 10,
"random_state": 1,
}

#create list to hold SSE values for each k
sse = []
for k in range(1, 11):
    kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
    kmeans.fit(principalupstreamwinnumscale)
    sse.append(kmeans.inertia_)

#visualize results
plt.plot(range(1, 11), sse)
plt.xticks(range(1, 11))
plt.xlabel("Number of Clusters")
plt.ylabel("SSE")
plt.show()


# In[56]:


principalupstreamwinnumscaledf=pd.DataFrame(principalupstreamwinnumscale)


# In[66]:


#append cluster assingments to original DataFrame\
principalupstreamwinnumscaledf=pd.DataFrame(principalupstreamwinnumscale)
principalupstreamwinnumscaledf['cluster'] = kmeans.labels_

#view updated DataFrame
print(principalupstreamwinnumscaledf)


# # Kmean PCA coloring Upstream

# In[29]:


for key in translater.keys():

    principalupstreamwinnumscaledf[key] = sumTRUP[key]
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.scatterplot(data=principalupstreamwinnumscaledf, x=principalupstreamwinnumscaledf['0'], y=principalupstreamwinnumscaledf['1'], hue=key)
    plt.title(key)
    plt.show()


# In[67]:


sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principalupstreamwinnumscaledf, x=principalupstreamwinnumscaledf[0], y=principalupstreamwinnumscaledf[1], hue='cluster')
plt.title("pca with ch38, kmean")
plt.show()


# In[20]:


## Save as TSV // read as CSV


# In[39]:


#principalupstreamwinnumscaledf.to_csv("principalupstreamwinnumscaledf", index= None)
principalupstreamwinnumscaledf=pd.read_csv('principalupstreamwinnumscaledf')
#upscaledff.to_csv("upscaledff", index= None)
upscaledff=pd.read_csv('upscaledff')


# # Kmean Upstream TSNE

# In[54]:


tsneupscaledf=StandardScaler().fit_transform(X_2u)
#initialize kmeans parameters
kmeans_kwargs = {
"init": "random",
"n_init": 10,
"random_state": 1,
}

#create list to hold SSE values for each k
sse = []
for k in range(1, 11):
    kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
    kmeans.fit(tsneupscaledf)
    sse.append(kmeans.inertia_)

#visualize results
plt.plot(range(1, 11), sse)
plt.xticks(range(1, 11))
plt.xlabel("Number of Clusters")
plt.ylabel("SSE")
plt.show()


# In[60]:


tsneupscaledf['cluster'] = kmeans.labels_


# In[59]:


tsneupscaledf=pd.DataFrame(tsneupscale)


# In[61]:


tsneupscaledf


# In[69]:


tsneupscaledf['cluster'].unique()


# In[71]:


tsneupscaledf['cluster'].value_counts()


# # Kmean TSNE coloring Upstream 

# In[64]:


sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=tsneupscaledf, x=tsneupscaledf[0], y=tsneupscaledf[1], hue='cluster')
plt.title("tsne, kmean")
plt.show()


# In[101]:


sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=tsneupscaledf[tsneupscaledf['cluster']== 2], x=tsneupscaledf[tsneupscaledf['cluster']== 2][0], y=tsneupscaledf[tsneupscaledf['cluster']== 2][1], hue='cluster')
plt.title("tsne, kmean")
plt.show()


# In[102]:


sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=tsneupscaledf[tsneupscaledf['cluster']== 8], x=tsneupscaledf[tsneupscaledf['cluster']== 8][0], y=tsneupscaledf[tsneupscaledf['cluster']== 8][1], hue='cluster')
plt.title("tsne, kmean")
plt.show()


# In[103]:


sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=tsneupscaledf[tsneupscaledf['cluster']== 0], x=tsneupscaledf[tsneupscaledf['cluster']== 0][0], y=tsneupscaledf[tsneupscaledf['cluster']== 0][1], hue='cluster')
plt.title("tsne, kmean")
plt.show()


# In[104]:


sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=tsneupscaledf[tsneupscaledf['cluster']== 9], x=tsneupscaledf[tsneupscaledf['cluster']== 9][0], y=tsneupscaledf[tsneupscaledf['cluster']== 9][1], hue='cluster')
plt.title("tsne, kmean")
plt.show()


# # UMAP Upstream

# In[111]:


# on UMAP
dfsumTRpltupscale=StandardScaler().fit_transform(dfsumTRpltup)
#initialize kmeans parameters
kmeans_kwargs = {
"init": "random",
"n_init": 10,
"random_state": 1,
}

#create list to hold SSE values for each k
sse = []
for k in range(1, 11):
    kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
    kmeans.fit(dfsumTRpltupscale)
    sse.append(kmeans.inertia_)

#visualize results
plt.plot(range(1, 11), sse)
plt.xticks(range(1, 11))
plt.xlabel("Number of Clusters")
plt.ylabel("SSE")
plt.show()


# In[ ]:


dfsumTRpltupscalef=pd.DataFrame(dfsumTRpltupscale)


# # Kmean UMAP coloring Upstream

# In[ ]:


for key in translater.keys():

    dfsumTRpltupscalef[key] = sumTRUP[key]
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.scatterplot(data= dfsumTRpltupscalef, x="a", y="b", hue=key)
    plt.title(key)
    plt.show()


# In[ ]:


dfsumTRpltupscalef['cluster'] = ch38drop['cluster']
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=dfsumTRpltupscalef, x="a", y="b", hue='cluster')
plt.title("UMAP with ch38, kmean")
plt.show()


# # Downstream Kmeans

# In[80]:


downstreamdropwinscale=StandardScaler().fit_transform(Downstreamdropwin)
#initialize kmeans parameters
kmeans_kwargs = {
"init": "random",
"n_init": 10,
"random_state": 1,
}

#create list to hold SSE values for each k
sse = []
for k in range(1, 11):
    kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
    kmeans.fit(downstreamdropwinscale)
    sse.append(kmeans.inertia_)

#visualize results
plt.plot(range(1, 11), sse)
plt.xticks(range(1, 11))
plt.xlabel("Number of Clusters")
plt.ylabel("SSE")
plt.show()


# In[85]:


#instantiate the k-means class, using optimal number of clusters
kmeans = KMeans(init="random", n_clusters=3, n_init=10, random_state=1)

#fit k-means algorithm to data
kmeans.fit(downstreamdropwinscale)

#view cluster assignments for each observation
kmeans.labels_

downscaledff=pd.DataFrame(downstreamdropwinscale)
#append cluster assingments to original DataFrame
downscaledff['cluster'] = kmeans.labels_

downscaledff.to_csv("downscaledff", index= None)
downscaledff
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=downscaledff, x=downscaledff[0], y=downscaledff[1], hue='cluster')
plt.title("tsne, kmean")
plt.show()


# # Downstream PCA kmeans 

# In[86]:


# on PCA
principaldownstreamwinnumscale=StandardScaler().fit_transform(principaldownstreamwinnum)
#initialize kmeans parameters
kmeans_kwargs = {
"init": "random",
"n_init": 10,
"random_state": 1,
}

#create list to hold SSE values for each k
sse = []
for k in range(1, 11):
    kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
    kmeans.fit(principaldownstreamwinnumscale)
    sse.append(kmeans.inertia_)

#visualize results
plt.plot(range(1, 11), sse)
plt.xticks(range(1, 11))
plt.xlabel("Number of Clusters")
plt.ylabel("SSE")
plt.show()


# In[89]:


#append cluster assingments to original DataFrame\
principaldownstreamwinnumscaledf=pd.DataFrame(principaldownstreamwinnumscale)
principaldownstreamwinnumscaledf['cluster'] = kmeans.labels_

#view updated DataFrame
print(principaldownstreamwinnumscaledf)


# In[90]:


sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=principaldownstreamwinnumscaledf, x=principaldownstreamwinnumscaledf[0], y=principaldownstreamwinnumscaledf[1], hue='cluster')
plt.title("pca, kmean")
plt.show()


# # DownStream Tsne Kmeans

# In[91]:


tsnedownscaledf=StandardScaler().fit_transform(X_2d)
#initialize kmeans parameters
kmeans_kwargs = {
"init": "random",
"n_init": 10,
"random_state": 1,
}

#create list to hold SSE values for each k
sse = []
for k in range(1, 11):
    kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
    kmeans.fit(tsnedownscaledf)
    sse.append(kmeans.inertia_)

#visualize results
plt.plot(range(1, 11), sse)
plt.xticks(range(1, 11))
plt.xlabel("Number of Clusters")
plt.ylabel("SSE")
plt.show()


# In[94]:


tsnedownscaledf=pd.DataFrame(tsnedownscaledf)


# In[95]:


tsnedownscaledf['cluster'] = kmeans.labels_


# In[96]:


sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.scatterplot(data=tsnedownscaledf, x=tsnedownscaledf[0], y=tsnedownscaledf[1], hue='cluster')
plt.title("tsne, kmean")
plt.show()


# In[ ]:


# on UMAP
dfsumTRpltdownscale=StandardScaler().fit_transform(dfsumTRpltdown)
#initialize kmeans parameters
kmeans_kwargs = {
"init": "random",
"n_init": 10,
"random_state": 1,
}

#create list to hold SSE values for each k
sse = []
for k in range(1, 11):
    kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
    kmeans.fit(dfsumTRpltdownscale)
    sse.append(kmeans.inertia_)

#visualize results
plt.plot(range(1, 11), sse)
plt.xticks(range(1, 11))
plt.xlabel("Number of Clusters")
plt.ylabel("SSE")
plt.show()


# # Louvain Community detection Algorithmn

# In[ ]:


import networkx as nx
import community as community_louvain
import matplotlib.pyplot as plt
import matplotlib.cm as cm
# define the graph
edge = [(1,2),(1,3),(1,4),(1,5),(1,6),(2,7),(2,8),(2,9)]
G = nx.Graph()
G.add_edges_from(edge)
# retrun partition as a dict
partition = community_louvain.best_partition(G)
# visualization
pos = nx.spring_layout(G)
cmap = cm.get_cmap('viridis', max(partition.values()) + 1)
nx.draw_networkx_nodes(G, pos, partition.keys(), node_size=100,cmap=cmap, node_color=list(partition.values()))
nx.draw_networkx_edges(G, pos, alpha=0.5)
plt.show()

