#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn import tree
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import make_pipeline
import seaborn as sns
import sklearn
#print(sklearn.__version__)
import graphviz
import pydotplus


# In[2]:


# read from master reduced excel file
df_all = pd.read_excel('Data/Degradation database v16_reduced.xlsx', index_col=0, header=0)
pd.set_option('display.max_rows', 1000)
df_all.columns


# In[3]:


# add colors column as heatmap for degredation
df_all['color_cat5'] = df_all['bio_rank_num5'].astype('category')
print(df_all['color_cat5'].cat.categories)
df_all['color_cat5'].cat.categories = ['red','orange','yellow','green','blue']

df_all['color_cat3'] = df_all['bio_rank_num3'].astype('category')
print(df_all['color_cat3'].cat.categories)
df_all['color_cat3'].cat.categories = [(0.267004, 0.004874, 0.329415),(0.128729, 0.563265, 0.551229),
                                       (0.993248, 0.906157, 0.143936)]
df_all['color_cat3_b'] = df_all['bio_rank_num3_b'].astype('category')
print(df_all['color_cat3_b'].cat.categories)
df_all['color_cat3_b'].cat.categories = ['orange','yellow','green']


# In[4]:


# sort by polymer type, a - abiotic, b - biotic
polyamidea = [84]
polyamideb = [78,83,85,86]

polycarbonatea = []
polycarbonateb = [79,80,81,82]

polyester_brancheda = [14,15,16] 
polyester_branchedb = [9,10,13,21,22,23,27,29,32,33,41,42,43,44] 
# ? 28b, 30b, 31b

polyester_cyclica = [17,18]
polyester_cyclicb = [11,75,76,77]
# ? 39a(PU), 74b

polyester_lineara = [19,25,26,48,50,52,54,56,58,60,62]
polyester_linearb = [0,1,2,3,4,5,6,7,8,12,20,24,34,35,36,37,38,45,46,47,49,51,53,55,57,59,61,63,64,65,66,67,68,69,70,71,72]
# ? 40b

polyethera = []
polyetherb = [90,91,92,93]

polyola = []
polyolb = [96,97]

polyolefina= []
polyolefinb= [88,89]

polyurethanea = [39]
polyurethaneb = [73]

vinylpyrrolidonea = []
vinylpyrrolidoneb = [94,95]

polystyrenea = []
polystyreneb = [99]

pvca = []
pvcb = [100]


# In[5]:


#try 3 cat equal classification trees bio - Tg Mn
rows = (polyamideb + polycarbonateb + polyester_branchedb + polyester_cyclicb + polyester_linearb + polyetherb +
       polyolb + polyolefinb + polyurethaneb + vinylpyrrolidoneb + polystyreneb + pvcb)
cols = ['Mn (kg/mol)','Tg (°C)','bio_rank_num3']
df_train = df_all.loc[rows,cols]
# remove NANs
df_train = df_train.drop([70,80,89,100])
#display(df_train)

# training set
x_train = df_train[cols[:-1]]
y_train = df_train['bio_rank_num3']

#train
treeb = tree.DecisionTreeClassifier(criterion='gini',max_depth=2,random_state=1)
fitb = treeb.fit(x_train,y_train)


# In[6]:


tree_g = tree.export_graphviz(fitb, out_file='tree4.dot', feature_names=['Mn (kg/mol)','Tg (°C)']
                              ,class_names=['slow','medium','fast'],impurity=False,rounded=True)
#graph = pydotplus.graph_from_dot_data(tree_g)
#graph.write_png('tree4.png')
#graph = graphviz.Source(tree_g)
#graph.view()


# In[7]:


# model evaluation
print('Number of samples: ' + str(df_train.index.size))
print('Accuracy score - full training set: %.3f' % fitb.score(x_train,y_train))
y_pred = treeb.predict(x_train)
conf_mat = confusion_matrix(y_true=y_train, y_pred=y_pred)
#print(conf_mat)
fig, ax = plt.subplots(figsize=(3.5,3.5))
ax.matshow(conf_mat, cmap=plt.cm.Blues, alpha=.5)
ax.set_xticklabels(['','Slow', 'Medium' ,'Fast'])
ax.set_yticklabels(['','Slow', 'Medium' ,'Fast'])
plt.xlabel('Predicted Label')
plt.ylabel('True Label')
for i in range(conf_mat.shape[0]):
    for j in range(conf_mat.shape[1]):
        ax.text(x=j,y=i,s=conf_mat[i,j],va='center',ha='center')
plt.tight_layout()

# cross-fold validation
treep = tree.DecisionTreeClassifier(criterion='gini',max_depth=2,random_state=1)
scores = cross_val_score(estimator=treep, X=x_train, y=y_train, cv=10)
#print(scores)
print('CV Accuracy Score: %.3f +/- %.3f' % (np.mean(scores),np.std(scores)))


# In[8]:


y_pred = treeb.predict(x_train)
correct = (y_pred != y_train.values)
df_incorrect = df_train[correct]
print(df_incorrect)


# In[9]:


# plot decision regions
plt.figure(figsize=(5.5,3.5),dpi=300)
rows = (polyamideb + polycarbonateb + polyester_branchedb + polyester_cyclicb + polyester_linearb + polyetherb +
       polyolb + polyolefinb + polyurethaneb + vinylpyrrolidoneb + polystyreneb + pvcb)
cols = ['Mn (kg/mol)','Tg (°C)','bio_rank_num3','color_cat3']
df_train = df_all.loc[rows,cols]
#display(df_train)
# linear polyester
Y = df_train.loc[polyester_linearb,'Mn (kg/mol)'].values
X = df_train.loc[polyester_linearb,'Tg (°C)'].values
plt_colors = df_train.loc[polyester_linearb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyester_linearb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='o')
# branched polyester
Y = df_train.loc[polyester_branchedb,'Mn (kg/mol)'].values
X = df_train.loc[polyester_branchedb,'Tg (°C)'].values
plt_colors = df_train.loc[polyester_branchedb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyester_branchedb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black', marker='p')
# cyclic polyester
Y = df_train.loc[polyester_cyclicb,'Mn (kg/mol)'].values
X = df_train.loc[polyester_cyclicb,'Tg (°C)'].values
plt_colors = df_train.loc[polyester_cyclicb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyester_cyclicb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='h')
# polyamide
Y = df_train.loc[polyamideb,'Mn (kg/mol)'].values
X = df_train.loc[polyamideb,'Tg (°C)'].values
plt_colors = df_train.loc[polyamideb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyamideb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black', marker='s')
# polycarbonate
Y = df_train.loc[polycarbonateb,'Mn (kg/mol)'].values
X = df_train.loc[polycarbonateb,'Tg (°C)'].values
plt_colors = df_train.loc[polycarbonateb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polycarbonateb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='*')
# polyether
Y = df_train.loc[polyetherb,'Mn (kg/mol)'].values
X = df_train.loc[polyetherb,'Tg (°C)'].values
plt_colors = df_train.loc[polyetherb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyetherb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='v')
# polyol
Y = df_train.loc[polyolb,'Mn (kg/mol)'].values
X = df_train.loc[polyolb,'Tg (°C)'].values
plt_colors = df_train.loc[polyolb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyolb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='<')
# polyolefin
Y = df_train.loc[polyolefinb,'Mn (kg/mol)'].values
X = df_train.loc[polyolefinb,'Tg (°C)'].values
plt_colors = df_train.loc[polyolefinb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyolefinb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='>')
# polyurethane
Y = df_train.loc[polyurethaneb,'Mn (kg/mol)'].values
X = df_train.loc[polyurethaneb,'Tg (°C)'].values
plt_colors = df_train.loc[polyurethaneb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyurethaneb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='^')
# vinylpyrrolidone
Y = df_train.loc[vinylpyrrolidoneb,'Mn (kg/mol)'].values
X = df_train.loc[vinylpyrrolidoneb,'Tg (°C)'].values
plt_colors = df_train.loc[vinylpyrrolidoneb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[vinylpyrrolidoneb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='D')
# polystyrene
Y = df_train.loc[polystyreneb,'Mn (kg/mol)'].values
X = df_train.loc[polystyreneb,'Tg (°C)'].values
plt_colors = df_train.loc[polystyreneb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polystyreneb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='P')
# pvc
Y = df_train.loc[pvcb,'Mn (kg/mol)'].values
X = df_train.loc[pvcb,'Tg (°C)'].values
plt_colors = df_train.loc[pvcb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[pvcb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='d')

# plot errors
#Y = df_incorrect['Mn (kg/mol)']
#X = df_incorrect['Tg (°C)']
#plt.scatter(X,Y,c='r',linewidths=1,marker='x')

# plot decision map
xx, yy, = np.meshgrid(np.arange(-110,201,1),np.arange(0,361.0,1.0))
z = fitb.predict(np.array([yy.ravel(),xx.ravel()]).T)
z = z.reshape(xx.shape)
#print(z)
cmap = mpl.colors.ListedColormap(((0.267004, 0.004874, 0.329415),(0.128729, 0.563265, 0.551229),
                                       (0.993248, 0.906157, 0.143936)))
plt.contourf(xx,yy,z,alpha=.2,cmap=cmap)

plt.ylabel('Mn (kg/mol)')
plt.xlabel('Tg (°C)')
plt.ylim(0,360.0)
plt.xlim(-110,200)

plt.tight_layout()
plt.savefig(fname='Fig4_122319.png',dpi=300,transparent=True,pad_inches=0.0)
#plt.show()


# In[27]:


# plot decision regions
plt.figure(figsize=(5.5,3.5),dpi=300)
rows = (polyamideb + polycarbonateb + polyester_branchedb + polyester_cyclicb + polyester_linearb + polyetherb +
       polyolb + polyolefinb + polyurethaneb + vinylpyrrolidoneb + polystyreneb + pvcb)
cols = ['Mn (kg/mol)','Tg (°C)','bio_rank_num3','color_cat3']
df_train = df_all.loc[rows,cols]
#display(df_train)
# linear polyester
Y = df_train.loc[polyester_linearb,'Mn (kg/mol)'].values
X = df_train.loc[polyester_linearb,'Tg (°C)'].values
plt_colors = df_train.loc[polyester_linearb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyester_linearb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='o')
# branched polyester
Y = df_train.loc[polyester_branchedb,'Mn (kg/mol)'].values
X = df_train.loc[polyester_branchedb,'Tg (°C)'].values
plt_colors = df_train.loc[polyester_branchedb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyester_branchedb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black', marker='p')
# cyclic polyester
Y = df_train.loc[polyester_cyclicb,'Mn (kg/mol)'].values
X = df_train.loc[polyester_cyclicb,'Tg (°C)'].values
plt_colors = df_train.loc[polyester_cyclicb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyester_cyclicb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='h')
# polyamide
Y = df_train.loc[polyamideb,'Mn (kg/mol)'].values
X = df_train.loc[polyamideb,'Tg (°C)'].values
plt_colors = df_train.loc[polyamideb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyamideb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black', marker='s')
# polycarbonate
Y = df_train.loc[polycarbonateb,'Mn (kg/mol)'].values
X = df_train.loc[polycarbonateb,'Tg (°C)'].values
plt_colors = df_train.loc[polycarbonateb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polycarbonateb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='*')
# polyether
Y = df_train.loc[polyetherb,'Mn (kg/mol)'].values
X = df_train.loc[polyetherb,'Tg (°C)'].values
plt_colors = df_train.loc[polyetherb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyetherb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='v')
# polyol
Y = df_train.loc[polyolb,'Mn (kg/mol)'].values
X = df_train.loc[polyolb,'Tg (°C)'].values
plt_colors = df_train.loc[polyolb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyolb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='<')
# polyolefin
Y = df_train.loc[polyolefinb,'Mn (kg/mol)'].values
X = df_train.loc[polyolefinb,'Tg (°C)'].values
plt_colors = df_train.loc[polyolefinb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyolefinb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='>')
# polyurethane
Y = df_train.loc[polyurethaneb,'Mn (kg/mol)'].values
X = df_train.loc[polyurethaneb,'Tg (°C)'].values
plt_colors = df_train.loc[polyurethaneb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyurethaneb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='^')
# vinylpyrrolidone
Y = df_train.loc[vinylpyrrolidoneb,'Mn (kg/mol)'].values
X = df_train.loc[vinylpyrrolidoneb,'Tg (°C)'].values
plt_colors = df_train.loc[vinylpyrrolidoneb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[vinylpyrrolidoneb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='D')
# polystyrene
Y = df_train.loc[polystyreneb,'Mn (kg/mol)'].values
X = df_train.loc[polystyreneb,'Tg (°C)'].values
plt_colors = df_train.loc[polystyreneb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polystyreneb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='P')
# pvc
Y = df_train.loc[pvcb,'Mn (kg/mol)'].values
X = df_train.loc[pvcb,'Tg (°C)'].values
plt_colors = df_train.loc[pvcb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[pvcb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='d')

# plot errors
Y = df_incorrect['Mn (kg/mol)']
X = df_incorrect['Tg (°C)']
plt.scatter(X,Y,c='black',linewidths=1,marker='x')

# plot decision map
xx, yy, = np.meshgrid(np.arange(-110,201,1),np.arange(0,361.0,1.0))
z = fitb.predict(np.array([yy.ravel(),xx.ravel()]).T)
z = z.reshape(xx.shape)
#print(z)
cmap = mpl.colors.ListedColormap(((0.267004, 0.004874, 0.329415),(0.128729, 0.563265, 0.551229),
                                       (0.993248, 0.906157, 0.143936)))
plt.contourf(xx,yy,z,alpha=.2,cmap=cmap)

plt.ylabel('$\mathit{M}\mathdefault{_{n}}$ $\mathdefault{(kg}$ $\mathdefault{mol^{-1})}$')
plt.xlabel('$\mathit{T}\mathdefault{_{g}}$ $\mathdefault{(°C)}$')
plt.ylim(0,360.0)
plt.xlim(-110,200)

plt.tight_layout()
#plt.savefig(fname='Fig4_122319.png',dpi=300,transparent=True,pad_inches=0.0)
plt.show()


# In[13]:


#try 3 cat equal classification trees bio - all colums - no cryst (corr w/ enth)
rows = (polyamideb + polycarbonateb + polyester_branchedb + polyester_cyclicb + polyester_linearb + polyetherb +
       polyolb + polyolefinb + polyurethaneb + vinylpyrrolidoneb + polystyreneb + pvcb)
cols = ['LogP/SA (Å-2)','Mn (kg/mol)','Tg (°C)','enthalpy (J/g)','bio_rank_num3']
df_train = df_all.loc[rows,cols]
# remove NANs
df_train = df_train.drop([10,11,12,36,38,44,70,78,80,89,94,95,99,100])
#display(df_train)

# training set
x_train = df_train[cols[:-1]]
y_train = df_train['bio_rank_num3']

#train
treeb = tree.DecisionTreeClassifier(criterion='gini',max_depth=3,random_state=1)
fitb = treeb.fit(x_train,y_train)


# In[14]:


tree_g = tree.export_graphviz(fitb, out_file='tree5.dot', feature_names=['LogP/SA (Å-2)','Mn (kg/mol)','Tg (°C)','enthalpy (J/g)']
                              ,class_names=['slow','medium','fast'],impurity=False,rounded=True)
#graph = graphviz.Source(tree_g)
#graph.view()


# In[15]:


# model evaluation
print('Number of samples: ' + str(df_train.index.size))
print('Accuracy score - full training set: %.3f' % fitb.score(x_train,y_train))
y_pred = treeb.predict(x_train)
conf_mat = confusion_matrix(y_true=y_train, y_pred=y_pred)
#print(conf_mat)
fig, ax = plt.subplots(figsize=(3.5,3.5))
ax.matshow(conf_mat, cmap=plt.cm.Blues, alpha=.5)
ax.set_xticklabels(['','Slow', 'Medium' ,'Fast'])
ax.set_yticklabels(['','Slow', 'Medium' ,'Fast'])
plt.xlabel('Predicted Label')
plt.ylabel('True Label')
for i in range(conf_mat.shape[0]):
    for j in range(conf_mat.shape[1]):
        ax.text(x=j,y=i,s=conf_mat[i,j],va='center',ha='center')
plt.tight_layout()

# cross-fold validation
treep = tree.DecisionTreeClassifier(criterion='gini',max_depth=2,random_state=1)
scores = cross_val_score(estimator=treep, X=x_train, y=y_train, cv=10)
#print(scores)
print('CV Accuracy Score: %.3f +/- %.3f' % (np.mean(scores),np.std(scores)))


# In[16]:


y_pred = treeb.predict(x_train)
correct = (y_pred != y_train.values)
df_incorrect = df_train[correct]
print(df_incorrect)


# In[17]:


#try 3 cat equal regression bio - LogPSA and Tg
rows = (polyamideb + polycarbonateb + polyester_branchedb + polyester_cyclicb + polyester_linearb + polyetherb +
       polyolb + polyolefinb + polyurethaneb + vinylpyrrolidoneb + polystyreneb + pvcb)
cols = ['LogP/SA (Å-2)','Tg (°C)','bio_rank_num3']
df_train = df_all.loc[rows,cols]
# remove NANs
df_train = df_train.drop([70])
#display(df_train)

# training set
x_train = df_train[cols[:-1]]
y_train = df_train['bio_rank_num3'].values

#standardization
std_sc = StandardScaler()
x_train_std = std_sc.fit_transform(x_train)
#display(x_train_std)

#train
svm = SVC(kernel='rbf', random_state=1, gamma=.2, C=10.0)
fitb = svm.fit(x_train_std,y_train)


# In[21]:


# plot decision regions
plt.figure(figsize=(5.75,3.5),dpi=300)
rows = (polyamideb + polycarbonateb + polyester_branchedb + polyester_cyclicb + polyester_linearb + polyetherb +
       polyolb + polyolefinb + polyurethaneb + vinylpyrrolidoneb + polystyreneb + pvcb)
cols = ['LogP/SA (Å-2)','Tg (°C)','bio_rank_num3','color_cat3']
df_train = df_all.loc[rows,cols]
#display(df_train)
# linear polyester
Y = df_train.loc[polyester_linearb,'LogP/SA (Å-2)'].values
X = df_train.loc[polyester_linearb,'Tg (°C)'].values
plt_colors = df_train.loc[polyester_linearb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyester_linearb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='o')
# branched polyester
Y = df_train.loc[polyester_branchedb,'LogP/SA (Å-2)'].values
X = df_train.loc[polyester_branchedb,'Tg (°C)'].values
plt_colors = df_train.loc[polyester_branchedb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyester_branchedb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black', marker='p')
# cyclic polyester
Y = df_train.loc[polyester_cyclicb,'LogP/SA (Å-2)'].values
X = df_train.loc[polyester_cyclicb,'Tg (°C)'].values
plt_colors = df_train.loc[polyester_cyclicb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyester_cyclicb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='h')
# polyamide
Y = df_train.loc[polyamideb,'LogP/SA (Å-2)'].values
X = df_train.loc[polyamideb,'Tg (°C)'].values
plt_colors = df_train.loc[polyamideb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyamideb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black', marker='s')
# polycarbonate
Y = df_train.loc[polycarbonateb,'LogP/SA (Å-2)'].values
X = df_train.loc[polycarbonateb,'Tg (°C)'].values
plt_colors = df_train.loc[polycarbonateb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polycarbonateb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='*')
# polyether
Y = df_train.loc[polyetherb,'LogP/SA (Å-2)'].values
X = df_train.loc[polyetherb,'Tg (°C)'].values
plt_colors = df_train.loc[polyetherb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyetherb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='v')
# polyol
Y = df_train.loc[polyolb,'LogP/SA (Å-2)'].values
X = df_train.loc[polyolb,'Tg (°C)'].values
plt_colors = df_train.loc[polyolb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyolb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='<')
# polyolefin
Y = df_train.loc[polyolefinb,'LogP/SA (Å-2)'].values
X = df_train.loc[polyolefinb,'Tg (°C)'].values
plt_colors = df_train.loc[polyolefinb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyolefinb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='>')
# polyurethane
Y = df_train.loc[polyurethaneb,'LogP/SA (Å-2)'].values
X = df_train.loc[polyurethaneb,'Tg (°C)'].values
plt_colors = df_train.loc[polyurethaneb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polyurethaneb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='^')
# vinylpyrrolidone
Y = df_train.loc[vinylpyrrolidoneb,'LogP/SA (Å-2)'].values
X = df_train.loc[vinylpyrrolidoneb,'Tg (°C)'].values
plt_colors = df_train.loc[vinylpyrrolidoneb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[vinylpyrrolidoneb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='D')
# polystyrene
Y = df_train.loc[polystyreneb,'LogP/SA (Å-2)'].values
X = df_train.loc[polystyreneb,'Tg (°C)'].values
plt_colors = df_train.loc[polystyreneb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[polystyreneb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='P')
# pvc
Y = df_train.loc[pvcb,'LogP/SA (Å-2)'].values
X = df_train.loc[pvcb,'Tg (°C)'].values
plt_colors = df_train.loc[pvcb,'color_cat3'].tolist()
plt_areas = (np.array(df_train.loc[pvcb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.7,linewidths=1,edgecolor='black',marker='d')

# plot decision map
# meshgrid is X(Tg),Y(LogPSA)
xx, yy = np.meshgrid(np.arange(-3.5,3.51,.01),np.arange(-3.5,3.51,.01))
#print([std_sc.transform(np.arange(-75,200,5)),std_sc.transform(np.arange(-.01,.025,.001))])
#xx_std, yy_std = np.meshgrid(std_sc.transform(np.arange(-75,200,5)),std_sc.transform(np.arange(-.01,.025,.001)))
z = fitb.predict(np.array([yy.ravel(),xx.ravel()]).T)
z = z.reshape(xx.shape)
test_df = pd.DataFrame()
test_df['LogPSA'] = np.arange(-3.5,3.51,.01)
test_df['Tg'] = np.arange(-3.5,3.51,.01)
test_df = std_sc.inverse_transform(test_df)
#print(test_df[:,0])
#print(z)
cmap = mpl.colors.ListedColormap(((0.267004, 0.004874, 0.329415),(0.128729, 0.563265, 0.551229),
                                       (0.993248, 0.906157, 0.143936)))
plt.contourf(test_df[:,1],test_df[:,0],z,alpha=.2,cmap=cmap)

# plot errors
Y = df_incorrect['LogP/SA (Å-2)']
X = df_incorrect['Tg (°C)']
plt.scatter(X,Y,c='white',linewidths=1,marker='x')

plt.ylabel('LogP/SA (Å-2)')
plt.xlabel('Tg (°C)')
plt.ylim(-.01,.025)
plt.xlim(-110,200)

plt.tight_layout()
plt.savefig(fname='Fig5_122319.png',dpi=300,transparent=True,pad_inches=0.0)
#plt.show()


# In[ ]:





# In[ ]:




