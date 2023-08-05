#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from matplotlib.colors import ListedColormap


# In[2]:


# read from master reduced excel file
df_all = pd.read_excel('Data/Degradation database v16_reduced.xlsx', index_col=0, header=0)
df_all.columns


# In[63]:


# add colors column as heatmap for degredation
df_all['color_cat5'] = df_all['bio_rank_num5'].astype('category')
print(df_all['color_cat5'].cat.categories)
df_all['color_cat5'].cat.categories = [(0.267004, 0.004874, 0.329415),(0.231674, 0.318106, 0.544834),
                                       (0.128729, 0.563265, 0.551229),(0.360741, 0.785964, 0.387814),
                                       (0.993248, 0.906157, 0.143936)]
df_all['color_cat3'] = df_all['bio_rank_num3'].astype('category')
print(df_all['color_cat3'].cat.categories)
df_all['color_cat3'].cat.categories = ['orange','yellow','green']
df_all['color_cat3_b'] = df_all['bio_rank_num3_b'].astype('category')
print(df_all['color_cat3_b'].cat.categories)
df_all['color_cat3_b'].cat.categories = ['orange','yellow','green']


# In[64]:


# Figure 2 polymers, a - abiotic, b - biotic
fig2a = [48,50,52,54,56,58,60,62,84]
fig2b = [47,49,51,53,55,57,59,61,83,85,86]
fig2 = fig2a + fig2b

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


# In[65]:


# Fig 2 - color map and area scaling
#color_map = plt.cm.get_cmap('jet_r')
color_map = plt.cm.get_cmap('viridis')
print(color_map.colors[0])
print(color_map.colors[63])
print(color_map.colors[127])
print(color_map.colors[191])
print(color_map.colors[255])

#print(color_map.size)

color_map = ListedColormap(color_map(np.linspace(0.0, 1.0, 256)))
scaler = MinMaxScaler(feature_range=(0.0,1.0))
print(np.min(np.array(df_all.loc[fig2,'wt. loss (mg/cm2 day)'].values)))
print(np.max(np.array(df_all.loc[fig2,'wt. loss (mg/cm2 day)'].values)))
deg_log_scale = np.log10(np.array(df_all.loc[fig2,'wt. loss (mg/cm2 day)'].values.reshape(-1,1)))
color_scale = scaler.fit_transform(deg_log_scale)
color_scale = color_scale[:,0]
plt_colors = [color_map(value) for value in color_scale]
#print(color_map.colors[0])
#print(color_map(19))
#print(color_map)

#print(plt_colors[127])
#print(plt_colors[191])
#print(plt_colors[255])
#print(plt_colors)
plt_sizes = deg_log_scale + 5.0
# note - slice the final arrays for the abio [:11] and bio [11:] data

# colorbar
fig, ax = plt.subplots(figsize=(.95, 3.0),dpi=300)
#fig.subplots_adjust(bottom=0.5)
#cmap = mpl.cm.cool
norm = mpl.colors.LogNorm(vmin=0.0001, vmax=100.0)
cb1 = mpl.colorbar.ColorbarBase(ax, cmap=color_map, norm=norm, 
                                orientation='vertical',extend='both')
#cb1.set_label('Degradation $(mg/cm^2/day)$\nSlow <----> Fast', fontsize=9)
plt.tick_params(labelsize=9)
#plt.clim(.0001,10)
plt.tight_layout()
plt.savefig(fname='Fig2scale_122319.png',dpi=300,transparent=True,pad_inches=0.0)
#fig.show()


# In[66]:


# Fig 2a - abiotic degradation, LogP/SA vs. Enthalpy of Melting
plt.figure(figsize=(3.25,2.5),dpi=300)
plt.scatter(df_all.loc[fig2a,'enthalpy (J/g)'].values,df_all.loc[fig2a,'LogP/SA (Å-2)'].values*1000,
            c=plt_colors[:9],s=(plt_sizes[:9]*1.4)**3.1,
            alpha=.95,linewidths=1,edgecolor='black')
plt.ylabel('LogP/SA $(Å^{-2})$  X $10^{-3}$', fontsize=9)
plt.xlabel('Enthalpy of Melting (J/g)', fontsize=9)
plt.tick_params(labelsize=9)
plt.ylim(0.0,13)
plt.xlim(39.0,84.0)
Xmod = [45.5,49.0,51,56,46.5,67.5,54,76,72]
Ymod = [0.0007,0.0034,0.0057,0.0077,0.0103,0.0108,0.0112,0.0105,0.005]
Ymod = [x * 1000 for x in Ymod]
#for i, txt in enumerate(df_all.loc[fig2a,'Name'].values):
#    plt.annotate(txt,(Xmod[i],(Ymod[i])),fontsize=9)
plt.tight_layout()
plt.savefig(fname='Fig2a_122319.png',dpi=300,transparent=True,pad_inches=0.0)
#plt.show()


# In[67]:


# Fig 2b - biotic degradation, LogP/SA vs. Enthalpy of Melting
plt.figure(figsize=(3.25,2.5),dpi=300)
plt.scatter(df_all.loc[fig2b,'enthalpy (J/g)'].values,df_all.loc[fig2b,'LogP/SA (Å-2)'].values*1000,
            c=plt_colors[9:],s=(plt_sizes[9:]*1.4)**3.1,
            alpha=.95,linewidths=1,edgecolor='black')
plt.ylabel('LogP/SA $(Å^{-2})$  X $10^{-3}$', fontsize=9)
plt.xlabel('Enthalpy of Melting (J/g)', fontsize=9)
plt.tick_params(labelsize=9)
plt.ylim(0.0,13)
plt.xlim(39.0,84.0)
Xmod = [45.5,50,40,44,45,68.5,53.5,77,74,75,64]
Ymod = [0.0004,0.0027,0.0062,0.008,0.01,0.011,0.0116,0.0096,0.0065,.0053,.0045]
Ymod = [x * 1000 for x in Ymod]
#for i, txt in enumerate(df_all.loc[fig2b,'Name'].values):
#    plt.annotate(txt,(Xmod[i],(Ymod[i])),fontsize=9)
plt.tight_layout()
plt.savefig(fname='Fig2b_122319.png',dpi=300,transparent=True,pad_inches=0.0)
#plt.show()


# In[68]:


# Fig 2c - abiotic degradation, LogP/SA vs. Crystallinity
plt.figure(figsize=(3.25,2.5),dpi=300)
plt.scatter(df_all.loc[fig2a,'% cryst'].values,df_all.loc[fig2a,'LogP/SA (Å-2)'].values*1000,
            c=plt_colors[:9],s=(plt_sizes[:9]*1.4)**3.1,
            alpha=.95,linewidths=1,edgecolor='black')
plt.ylabel('LogP/SA $(Å^{-2})$  X $10^{-3}$', fontsize=9)
plt.xlabel('Crystallinity (%)',fontsize=9)
plt.tick_params(labelsize=9)
plt.ylim(0.0,13)
plt.xlim(26.0,58.0)
Xmod = [28.2,26.5,31.0,39.5,34,47.0,45.5,51.5,35]
Ymod = [0.0007,0.004,0.0068,0.008,0.0104,0.0104,0.0116,0.0105,0.0035]
Ymod = [x * 1000 for x in Ymod]
#for i, txt in enumerate(df_all.loc[fig2a,'Name'].values):
#    plt.annotate(txt,(Xmod[i],(Ymod[i])),fontsize=9)
plt.tight_layout()
plt.savefig(fname='Fig2c_122319.png',dpi=300,transparent=True,pad_inches=0.0)
#plt.show()


# In[69]:


# Fig 2d - biotic degradation, LogP/SA vs. Cryst
plt.figure(figsize=(3.25,2.5),dpi=300)
plt.scatter(df_all.loc[fig2b,'% cryst'].values,df_all.loc[fig2b,'LogP/SA (Å-2)'].values*1000,
            c=plt_colors[9:],s=(plt_sizes[9:]*1.4)**3.1,
            alpha=.95,linewidths=1,edgecolor='black')
plt.ylabel('LogP/SA $(Å^{-2})$  X $10^{-3}$', fontsize=9)
plt.xlabel('Crystallinity (%)',fontsize=9)
plt.tick_params(labelsize=9)
plt.ylim(0.0,13)
plt.xlim(26.0,58.0)
Xmod = [27.5,26.5,30.0,41,32,47.5,38,55.,35.2,37,33]
Ymod = [0.0003,0.0044,0.0068,0.008,0.0104,0.0108,0.0117,0.0108,0.0034,.0043,.0025]
Ymod = [x * 1000 for x in Ymod]
#for i, txt in enumerate(df_all.loc[fig2b,'Name'].values):
#    plt.annotate(txt,(Xmod[i],(Ymod[i])),fontsize=9)
plt.tight_layout()
plt.savefig(fname='Fig2d_122319.png',dpi=300,transparent=True,pad_inches=0.0)
plt.show()


# In[74]:


### Fig 3a all polymers - abiotic - low Tg - based on 5-categorical degradation ranking
plt.figure(figsize=(3.42,2.5),dpi=300)
# linear polyester
Y = df_all.loc[polyester_lineara,'Mn (kg/mol)'].values
X = df_all.loc[polyester_lineara,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_lineara,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyester_lineara,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='o')
# branched polyester
Y = df_all.loc[polyester_brancheda,'Mn (kg/mol)'].values
X = df_all.loc[polyester_brancheda,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_brancheda,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyester_brancheda,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black', marker='p')
# cyclic polyester
Y = df_all.loc[polyester_cyclica,'Mn (kg/mol)'].values
X = df_all.loc[polyester_cyclica,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_cyclica,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyester_cyclica,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='h')
# polyamide
Y = df_all.loc[polyamidea,'Mn (kg/mol)'].values
X = df_all.loc[polyamidea,'Tg (°C)'].values
plt_colors = df_all.loc[polyamidea,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyamidea,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black', marker='s')
# polycarbonate
Y = df_all.loc[polycarbonatea,'Mn (kg/mol)'].values
X = df_all.loc[polycarbonatea,'Tg (°C)'].values
plt_colors = df_all.loc[polycarbonatea,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polycarbonatea,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='*')
# polyether
Y = df_all.loc[polyethera,'Mn (kg/mol)'].values
X = df_all.loc[polyethera,'Tg (°C)'].values
plt_colors = df_all.loc[polyethera,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyethera,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='v')
# polyol
Y = df_all.loc[polyola,'Mn (kg/mol)'].values
X = df_all.loc[polyola,'Tg (°C)'].values
plt_colors = df_all.loc[polyola,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyola,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='<')
# polyolefin
Y = df_all.loc[polyolefina,'Mn (kg/mol)'].values
X = df_all.loc[polyolefina,'Tg (°C)'].values
plt_colors = df_all.loc[polyolefina,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyolefina,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='>')
# polyurethane
Y = df_all.loc[polyurethanea,'Mn (kg/mol)'].values
X = df_all.loc[polyurethanea,'Tg (°C)'].values
plt_colors = df_all.loc[polyurethanea,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyurethanea,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='^')
# vinylpyrrolidone
Y = df_all.loc[vinylpyrrolidonea,'Mn (kg/mol)'].values
X = df_all.loc[vinylpyrrolidonea,'Tg (°C)'].values
plt_colors = df_all.loc[vinylpyrrolidonea,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[vinylpyrrolidonea,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='D')
# polystyrene
Y = df_all.loc[polystyrenea,'Mn (kg/mol)'].values
X = df_all.loc[polystyrenea,'Tg (°C)'].values
plt_colors = df_all.loc[polystyrenea,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polystyrenea,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='P')
# pvc
Y = df_all.loc[pvca,'Mn (kg/mol)'].values
X = df_all.loc[pvca,'Tg (°C)'].values
plt_colors = df_all.loc[pvca,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[pvca,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='d')

plt.ylabel('Mn (kg/mol)', fontsize=9)
plt.xlabel('Tg (°C)', fontsize=9)
plt.tick_params(labelsize=9)
plt.ylim(0,60)
plt.xlim(-70,0)
plt.tight_layout()
plt.savefig(fname='Fig3a_122319.png',dpi=300,transparent=True,pad_inches=0.0)
plt.show()


# In[75]:


# Fig 3c all polymers - abiotic - High Tg - based on 5-categorical degradation ranking
plt.figure(figsize=(3.5,2.5),dpi=300)
# linear polyester
Y = df_all.loc[polyester_lineara,'Mn (kg/mol)'].values
X = df_all.loc[polyester_lineara,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_lineara,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyester_lineara,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='o')
# branched polyester
Y = df_all.loc[polyester_brancheda,'Mn (kg/mol)'].values
X = df_all.loc[polyester_brancheda,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_brancheda,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyester_brancheda,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black', marker='p')
# cyclic polyester
Y = df_all.loc[polyester_cyclica,'Mn (kg/mol)'].values
X = df_all.loc[polyester_cyclica,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_cyclica,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyester_cyclica,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='h')
# polyamide
Y = df_all.loc[polyamidea,'Mn (kg/mol)'].values
X = df_all.loc[polyamidea,'Tg (°C)'].values
plt_colors = df_all.loc[polyamidea,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyamidea,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black', marker='s')
# polycarbonate
Y = df_all.loc[polycarbonatea,'Mn (kg/mol)'].values
X = df_all.loc[polycarbonatea,'Tg (°C)'].values
plt_colors = df_all.loc[polycarbonatea,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polycarbonatea,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='*')
# polyether
Y = df_all.loc[polyethera,'Mn (kg/mol)'].values
X = df_all.loc[polyethera,'Tg (°C)'].values
plt_colors = df_all.loc[polyethera,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyethera,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='v')
# polyol
Y = df_all.loc[polyola,'Mn (kg/mol)'].values
X = df_all.loc[polyola,'Tg (°C)'].values
plt_colors = df_all.loc[polyola,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyola,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='<')
# polyolefin
Y = df_all.loc[polyolefina,'Mn (kg/mol)'].values
X = df_all.loc[polyolefina,'Tg (°C)'].values
plt_colors = df_all.loc[polyolefina,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyolefina,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='>')
# polyurethane
Y = df_all.loc[polyurethanea,'Mn (kg/mol)'].values
X = df_all.loc[polyurethanea,'Tg (°C)'].values
plt_colors = df_all.loc[polyurethanea,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyurethanea,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='^')
# vinylpyrrolidone
Y = df_all.loc[vinylpyrrolidonea,'Mn (kg/mol)'].values
X = df_all.loc[vinylpyrrolidonea,'Tg (°C)'].values
plt_colors = df_all.loc[vinylpyrrolidonea,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[vinylpyrrolidonea,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='D')
# polystyrene
Y = df_all.loc[polystyrenea,'Mn (kg/mol)'].values
X = df_all.loc[polystyrenea,'Tg (°C)'].values
plt_colors = df_all.loc[polystyrenea,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polystyrenea,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='P')
# pvc
Y = df_all.loc[pvca,'Mn (kg/mol)'].values
X = df_all.loc[pvca,'Tg (°C)'].values
plt_colors = df_all.loc[pvca,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[pvca,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='d')

plt.ylabel('Mn (kg/mol)', fontsize=9)
plt.xlabel('Tg (°C)', fontsize=9)
plt.tick_params(labelsize=9)
plt.ylim(0,380)
plt.xlim(-10,180)
plt.tight_layout()
plt.savefig(fname='Fig3c_122319.png',dpi=300,transparent=True,pad_inches=0.0)
plt.show()


# In[76]:


# Fig 3b all polymers - biotic - low Tg - based on 5-categorical degradation ranking
plt.figure(figsize=(3.5,2.5),dpi=300)
# linear polyester
Y = df_all.loc[polyester_linearb,'Mn (kg/mol)'].values
X = df_all.loc[polyester_linearb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_linearb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyester_linearb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='o')
# branched polyester
Y = df_all.loc[polyester_branchedb,'Mn (kg/mol)'].values
X = df_all.loc[polyester_branchedb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_branchedb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyester_branchedb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black', marker='p')
# cyclic polyester
Y = df_all.loc[polyester_cyclicb,'Mn (kg/mol)'].values
X = df_all.loc[polyester_cyclicb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_cyclicb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyester_cyclicb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='h')
# polyamide
Y = df_all.loc[polyamideb,'Mn (kg/mol)'].values
X = df_all.loc[polyamideb,'Tg (°C)'].values
plt_colors = df_all.loc[polyamideb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyamideb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black', marker='s')
# polycarbonate
Y = df_all.loc[polycarbonateb,'Mn (kg/mol)'].values
X = df_all.loc[polycarbonateb,'Tg (°C)'].values
plt_colors = df_all.loc[polycarbonateb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polycarbonateb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='*')
# polyether
Y = df_all.loc[polyetherb,'Mn (kg/mol)'].values
X = df_all.loc[polyetherb,'Tg (°C)'].values
plt_colors = df_all.loc[polyetherb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyetherb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='v')
# polyol
Y = df_all.loc[polyolb,'Mn (kg/mol)'].values
X = df_all.loc[polyolb,'Tg (°C)'].values
plt_colors = df_all.loc[polyolb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyolb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='<')
# polyolefin
Y = df_all.loc[polyolefinb,'Mn (kg/mol)'].values
X = df_all.loc[polyolefinb,'Tg (°C)'].values
plt_colors = df_all.loc[polyolefinb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyolefinb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='>')
# polyurethane
Y = df_all.loc[polyurethaneb,'Mn (kg/mol)'].values
X = df_all.loc[polyurethaneb,'Tg (°C)'].values
plt_colors = df_all.loc[polyurethaneb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyurethaneb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='^')
# vinylpyrrolidone
Y = df_all.loc[vinylpyrrolidoneb,'Mn (kg/mol)'].values
X = df_all.loc[vinylpyrrolidoneb,'Tg (°C)'].values
plt_colors = df_all.loc[vinylpyrrolidoneb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[vinylpyrrolidoneb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='D')
# polystyrene
Y = df_all.loc[polystyreneb,'Mn (kg/mol)'].values
X = df_all.loc[polystyreneb,'Tg (°C)'].values
plt_colors = df_all.loc[polystyreneb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polystyreneb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='P')
# pvc
Y = df_all.loc[pvcb,'Mn (kg/mol)'].values
X = df_all.loc[pvcb,'Tg (°C)'].values
plt_colors = df_all.loc[pvcb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[pvcb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='d')

plt.ylabel('Mn (kg/mol)', fontsize=9)
plt.xlabel('Tg (°C)', fontsize=9)
plt.tick_params(labelsize=9)
plt.ylim(0,120)
plt.xlim(-70,0)
plt.tight_layout()
plt.savefig(fname='Fig3b_122319.png',dpi=300,transparent=True,pad_inches=0.0)
plt.show()


# In[77]:


# Fig 3d all polymers - biotic - High Tg - based on 5-categorical degradation ranking
plt.figure(figsize=(3.5,2.5),dpi=300)
# linear polyester
Y = df_all.loc[polyester_linearb,'Mn (kg/mol)'].values
X = df_all.loc[polyester_linearb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_linearb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyester_linearb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='o')
# branched polyester
Y = df_all.loc[polyester_branchedb,'Mn (kg/mol)'].values
X = df_all.loc[polyester_branchedb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_branchedb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyester_branchedb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black', marker='p')
# cyclic polyester
Y = df_all.loc[polyester_cyclicb,'Mn (kg/mol)'].values
X = df_all.loc[polyester_cyclicb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_cyclicb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyester_cyclicb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='h')
# polyamide
Y = df_all.loc[polyamideb,'Mn (kg/mol)'].values
X = df_all.loc[polyamideb,'Tg (°C)'].values
plt_colors = df_all.loc[polyamideb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyamideb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black', marker='s')
# polycarbonate
Y = df_all.loc[polycarbonateb,'Mn (kg/mol)'].values
X = df_all.loc[polycarbonateb,'Tg (°C)'].values
plt_colors = df_all.loc[polycarbonateb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polycarbonateb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='*')
# polyether
Y = df_all.loc[polyetherb,'Mn (kg/mol)'].values
X = df_all.loc[polyetherb,'Tg (°C)'].values
plt_colors = df_all.loc[polyetherb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyetherb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='v')
# polyol
Y = df_all.loc[polyolb,'Mn (kg/mol)'].values
X = df_all.loc[polyolb,'Tg (°C)'].values
plt_colors = df_all.loc[polyolb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyolb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='<')
# polyolefin
Y = df_all.loc[polyolefinb,'Mn (kg/mol)'].values
X = df_all.loc[polyolefinb,'Tg (°C)'].values
plt_colors = df_all.loc[polyolefinb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyolefinb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='>')
# polyurethane
Y = df_all.loc[polyurethaneb,'Mn (kg/mol)'].values
X = df_all.loc[polyurethaneb,'Tg (°C)'].values
plt_colors = df_all.loc[polyurethaneb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyurethaneb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='^')
# vinylpyrrolidone
Y = df_all.loc[vinylpyrrolidoneb,'Mn (kg/mol)'].values
X = df_all.loc[vinylpyrrolidoneb,'Tg (°C)'].values
plt_colors = df_all.loc[vinylpyrrolidoneb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[vinylpyrrolidoneb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='D')
# polystyrene
Y = df_all.loc[polystyreneb,'Mn (kg/mol)'].values
X = df_all.loc[polystyreneb,'Tg (°C)'].values
plt_colors = df_all.loc[polystyreneb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polystyreneb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='P')
# pvc
Y = df_all.loc[pvcb,'Mn (kg/mol)'].values
X = df_all.loc[pvcb,'Tg (°C)'].values
plt_colors = df_all.loc[pvcb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[pvcb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.8,linewidths=1,edgecolor='black',marker='d')

plt.ylabel('Mn (kg/mol)', fontsize=9)
plt.xlabel('Tg (°C)', fontsize=9)
plt.tick_params(labelsize=9)
plt.ylim(0,380)
plt.xlim(-10,180)
plt.tight_layout()
plt.savefig(fname='Fig3d_122319.png',dpi=300,transparent=True,pad_inches=0.0)
plt.show()


# In[30]:


# Fig 3e polyester branched and cyclic and other cyclics - abiotic - High Tg - based on 5-categorical degradation ranking
plt.figure(figsize=(3.5,2.5),dpi=100)

# branched polyester
Y = df_all.loc[polyester_brancheda,'Mn (kg/mol)'].values
X = df_all.loc[polyester_brancheda,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_brancheda,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyester_brancheda,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black', marker='p')
# cyclic polyester
Y = df_all.loc[polyester_cyclica,'Mn (kg/mol)'].values
X = df_all.loc[polyester_cyclica,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_cyclica,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyester_cyclica,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='h')
# polycarbonate
Y = df_all.loc[polycarbonatea,'Mn (kg/mol)'].values
X = df_all.loc[polycarbonatea,'Tg (°C)'].values
plt_colors = df_all.loc[polycarbonatea,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polycarbonatea,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='*')
# vinylpyrrolidone
Y = df_all.loc[vinylpyrrolidonea,'Mn (kg/mol)'].values
X = df_all.loc[vinylpyrrolidonea,'Tg (°C)'].values
plt_colors = df_all.loc[vinylpyrrolidonea,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[vinylpyrrolidonea,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='D')


plt.ylabel('Mn (kg/mol)', fontsize=9)
plt.xlabel('Tg (°C)', fontsize=9)
plt.tick_params(labelsize=9)
plt.ylim(0,360)
plt.xlim(-5,180)
plt.tight_layout()
plt.savefig(fname='Fig3e_070219.png',dpi=300,transparent=True,pad_inches=0.0)
plt.show()


# In[31]:


# Fig 3f polyester branched and cyclic and other cyclics - biotic - High Tg - based on 5-categorical degradation ranking
plt.figure(figsize=(3.5,2.5),dpi=100)

# branched polyester
Y = df_all.loc[polyester_branchedb,'Mn (kg/mol)'].values
X = df_all.loc[polyester_branchedb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_branchedb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyester_branchedb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black', marker='p')
# cyclic polyester
Y = df_all.loc[polyester_cyclicb,'Mn (kg/mol)'].values
X = df_all.loc[polyester_cyclicb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_cyclicb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polyester_cyclicb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='h')
# polycarbonate
Y = df_all.loc[polycarbonateb,'Mn (kg/mol)'].values
X = df_all.loc[polycarbonateb,'Tg (°C)'].values
plt_colors = df_all.loc[polycarbonateb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[polycarbonateb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='*')
# vinylpyrrolidone
Y = df_all.loc[vinylpyrrolidoneb,'Mn (kg/mol)'].values
X = df_all.loc[vinylpyrrolidoneb,'Tg (°C)'].values
plt_colors = df_all.loc[vinylpyrrolidoneb,'color_cat5'].tolist()
plt_areas = (np.array(df_all.loc[vinylpyrrolidoneb,'bio_rank_num5'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='D')

plt.ylabel('Mn (kg/mol)', fontsize=9)
plt.xlabel('Tg (°C)', fontsize=9)
plt.tick_params(labelsize=9)
plt.ylim(0,360)
plt.xlim(-5,180)
plt.tight_layout()
plt.savefig(fname='Fig3f_070219.png',dpi=300,transparent=True,pad_inches=0.0)
plt.show()


# In[5]:


# Legend
#elements = [([],[],marker='o',label='poly'),([],[],marker='*',label='carb')]
plt.figure(figsize=(3.5,4.5),dpi=300)
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='o',label='linear polyester')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='p',label='branched polyester')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='h',label='cyclic polyester')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='s',label='polyamide')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='*',label='polycarbonate')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='v',label='polyether')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='<',label='poly(vinyl alcohol)')
#plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='>',label='polyolefin')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='^',label='polyurethane')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='D',label='polyvinylpyrrolidone')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='P',label='polystyrene')
#plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='d',label='polyvinyl chloride')

plt.legend(fontsize=9,loc='center')
plt.tight_layout()
plt.savefig(fname='Fig3legend_122319.png',dpi=300,transparent=True,pad_inches=0.0)
plt.show()


# In[6]:


# Legend full
#elements = [([],[],marker='o',label='poly'),([],[],marker='*',label='carb')]
plt.figure(figsize=(3.5,4.5),dpi=300)
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='o',label='linear polyester')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='p',label='branched polyester')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='h',label='cyclic polyester')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='s',label='polyamide')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='*',label='polycarbonate')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='v',label='polyether')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='<',label='poly(vinyl alcohol)')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='>',label='polyolefin')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='^',label='polyurethane')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='D',label='polyvinylpyrrolidone')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='P',label='polystyrene')
plt.scatter([],[],c='white',s=50,alpha=1.0,linewidths=1,edgecolor='black',marker='d',label='polyvinyl chloride')

plt.legend(fontsize=9,loc='center')
plt.tight_layout()
plt.savefig(fname='Fig3legendfull_122319.png',dpi=300,transparent=True,pad_inches=0.0)
plt.show()


# In[26]:


# Fig 4 test - LogPSA versus Tg - bio - 3cat1 - all polymers
plt.figure(figsize=(7.5,5.5),dpi=100)
# linear polyester
Y = df_all.loc[polyester_linearb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyester_linearb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_linearb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polyester_linearb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='o')
# branched polyester
Y = df_all.loc[polyester_branchedb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyester_branchedb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_branchedb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polyester_branchedb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black', marker='p')
# cyclic polyester
Y = df_all.loc[polyester_cyclicb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyester_cyclicb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_cyclicb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polyester_cyclicb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='h')
# polyamide
Y = df_all.loc[polyamideb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyamideb,'Tg (°C)'].values
plt_colors = df_all.loc[polyamideb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polyamideb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black', marker='s')
# polycarbonate
Y = df_all.loc[polycarbonateb,'LogP/SA (Å-2)'].values
X = df_all.loc[polycarbonateb,'Tg (°C)'].values
plt_colors = df_all.loc[polycarbonateb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polycarbonateb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='*')
# polyether
Y = df_all.loc[polyetherb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyetherb,'Tg (°C)'].values
plt_colors = df_all.loc[polyetherb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polyetherb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='v')
# polyol
Y = df_all.loc[polyolb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyolb,'Tg (°C)'].values
plt_colors = df_all.loc[polyolb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polyolb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='<')
# polyolefin
Y = df_all.loc[polyolefinb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyolefinb,'Tg (°C)'].values
plt_colors = df_all.loc[polyolefinb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polyolefinb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='>')
# polyurethane
Y = df_all.loc[polyurethaneb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyurethaneb,'Tg (°C)'].values
plt_colors = df_all.loc[polyurethaneb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polyurethaneb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='^')
# vinylpyrrolidone
Y = df_all.loc[vinylpyrrolidoneb,'LogP/SA (Å-2)'].values
X = df_all.loc[vinylpyrrolidoneb,'Tg (°C)'].values
plt_colors = df_all.loc[vinylpyrrolidoneb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[vinylpyrrolidoneb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='D')
# polystyrene
Y = df_all.loc[polystyreneb,'LogP/SA (Å-2)'].values
X = df_all.loc[polystyreneb,'Tg (°C)'].values
plt_colors = df_all.loc[polystyreneb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polystyreneb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='P')
# pvc
Y = df_all.loc[pvcb,'LogP/SA (Å-2)'].values
X = df_all.loc[pvcb,'Tg (°C)'].values
plt_colors = df_all.loc[pvcb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[pvcb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='d')

plt.ylabel('LogP/SA',fontsize=9)
plt.xlabel('Tg (C)',fontsize=9)
plt.tick_params(labelsize=9)
plt.ylim(-.01,.025)
plt.xlim(-110,200)
plt.tight_layout()
#plt.savefig(fname='Fig3a_062419.png',dpi=300,transparent=True,pad_inches=0.0)
plt.show()


# In[24]:


# Fig 4 test - LogPSA versus Tg - bio - 3cat1 - all polymers
plt.figure(figsize=(7.5,5.5),dpi=100)
# linear polyester
Y = df_all.loc[polyester_linearb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyester_linearb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_linearb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polyester_linearb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='o')
# branched polyester
Y = df_all.loc[polyester_branchedb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyester_branchedb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_branchedb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polyester_branchedb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black', marker='p')
# cyclic polyester
Y = df_all.loc[polyester_cyclicb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyester_cyclicb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_cyclicb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polyester_cyclicb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='h')
# polyamide
Y = df_all.loc[polyamideb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyamideb,'Tg (°C)'].values
plt_colors = df_all.loc[polyamideb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polyamideb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black', marker='s')
# polycarbonate
Y = df_all.loc[polycarbonateb,'LogP/SA (Å-2)'].values
X = df_all.loc[polycarbonateb,'Tg (°C)'].values
plt_colors = df_all.loc[polycarbonateb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polycarbonateb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='*')
# polyether
Y = df_all.loc[polyetherb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyetherb,'Tg (°C)'].values
plt_colors = df_all.loc[polyetherb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polyetherb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='v')
# polyol
Y = df_all.loc[polyolb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyolb,'Tg (°C)'].values
plt_colors = df_all.loc[polyolb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polyolb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='<')
# polyolefin
Y = df_all.loc[polyolefinb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyolefinb,'Tg (°C)'].values
plt_colors = df_all.loc[polyolefinb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polyolefinb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='>')
# polyurethane
Y = df_all.loc[polyurethaneb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyurethaneb,'Tg (°C)'].values
plt_colors = df_all.loc[polyurethaneb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polyurethaneb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='^')
# vinylpyrrolidone
Y = df_all.loc[vinylpyrrolidoneb,'LogP/SA (Å-2)'].values
X = df_all.loc[vinylpyrrolidoneb,'Tg (°C)'].values
plt_colors = df_all.loc[vinylpyrrolidoneb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[vinylpyrrolidoneb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='D')
# polystyrene
Y = df_all.loc[polystyreneb,'LogP/SA (Å-2)'].values
X = df_all.loc[polystyreneb,'Tg (°C)'].values
plt_colors = df_all.loc[polystyreneb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[polystyreneb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='P')
# pvc
Y = df_all.loc[pvcb,'LogP/SA (Å-2)'].values
X = df_all.loc[pvcb,'Tg (°C)'].values
plt_colors = df_all.loc[pvcb,'color_cat3'].tolist()
plt_areas = (np.array(df_all.loc[pvcb,'bio_rank_num3'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='d')

plt.ylabel('LogP/SA',fontsize=9)
plt.xlabel('Tg (C)',fontsize=9)
plt.tick_params(labelsize=9)

all_plot = (polyester_linearb + polyester_branchedb + polyester_cyclicb + polyamideb + polycarbonateb + polyetherb +
            polyolb + polyolefinb + polyurethaneb + vinylpyrrolidoneb + polystyreneb + pvcb)
for i in np.nditer(df_all.loc[all_plot].index):
    try:
        Xmod2 = float(df_all.loc[int(i),'Tg (°C)'])
        if np.isnan(Xmod2):
            Xmod = -75
        else:
            Xmod = Xmod2
    except ValueError:
        print('here')
        Xmod = -75
    try:
        Ymod2 = float(df_all.loc[int(i),'LogP/SA (Å-2)'])
        if np.isnan(Ymod2):
            Ymod = -.01
        else:
            Ymod = Ymod2
    except ValueError:
        print('here2')
        Ymod = -.01
    plt.annotate(i,(Xmod-4.0,Ymod-.0004),fontsize=9)
    
plt.ylim(-.01,.025)
plt.xlim(-75,200)
plt.tight_layout()
#plt.savefig(fname='Fig3a_062419.png',dpi=300,transparent=True,pad_inches=0.0)
plt.show()


# In[25]:


# Fig 4 test - LogPSA versus Tg - bio - 3cat2 - all polymers
plt.figure(figsize=(7.5,5.5),dpi=100)
# linear polyester
Y = df_all.loc[polyester_linearb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyester_linearb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_linearb,'color_cat3_b'].tolist()
plt_areas = (np.array(df_all.loc[polyester_linearb,'bio_rank_num3_b'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='o')
# branched polyester
Y = df_all.loc[polyester_branchedb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyester_branchedb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_branchedb,'color_cat3_b'].tolist()
plt_areas = (np.array(df_all.loc[polyester_branchedb,'bio_rank_num3_b'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black', marker='p')
# cyclic polyester
Y = df_all.loc[polyester_cyclicb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyester_cyclicb,'Tg (°C)'].values
plt_colors = df_all.loc[polyester_cyclicb,'color_cat3_b'].tolist()
plt_areas = (np.array(df_all.loc[polyester_cyclicb,'bio_rank_num3_b'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='h')
# polyamide
Y = df_all.loc[polyamideb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyamideb,'Tg (°C)'].values
plt_colors = df_all.loc[polyamideb,'color_cat3_b'].tolist()
plt_areas = (np.array(df_all.loc[polyamideb,'bio_rank_num3_b'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black', marker='s')
# polycarbonate
Y = df_all.loc[polycarbonateb,'LogP/SA (Å-2)'].values
X = df_all.loc[polycarbonateb,'Tg (°C)'].values
plt_colors = df_all.loc[polycarbonateb,'color_cat3_b'].tolist()
plt_areas = (np.array(df_all.loc[polycarbonateb,'bio_rank_num3_b'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='*')
# polyether
Y = df_all.loc[polyetherb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyetherb,'Tg (°C)'].values
plt_colors = df_all.loc[polyetherb,'color_cat3_b'].tolist()
plt_areas = (np.array(df_all.loc[polyetherb,'bio_rank_num3_b'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='v')
# polyol
Y = df_all.loc[polyolb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyolb,'Tg (°C)'].values
plt_colors = df_all.loc[polyolb,'color_cat3_b'].tolist()
plt_areas = (np.array(df_all.loc[polyolb,'bio_rank_num3_b'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='<')
# polyolefin
Y = df_all.loc[polyolefinb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyolefinb,'Tg (°C)'].values
plt_colors = df_all.loc[polyolefinb,'color_cat3_b'].tolist()
plt_areas = (np.array(df_all.loc[polyolefinb,'bio_rank_num3_b'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='>')
# polyurethane
Y = df_all.loc[polyurethaneb,'LogP/SA (Å-2)'].values
X = df_all.loc[polyurethaneb,'Tg (°C)'].values
plt_colors = df_all.loc[polyurethaneb,'color_cat3_b'].tolist()
plt_areas = (np.array(df_all.loc[polyurethaneb,'bio_rank_num3_b'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='^')
# vinylpyrrolidone
Y = df_all.loc[vinylpyrrolidoneb,'LogP/SA (Å-2)'].values
X = df_all.loc[vinylpyrrolidoneb,'Tg (°C)'].values
plt_colors = df_all.loc[vinylpyrrolidoneb,'color_cat3_b'].tolist()
plt_areas = (np.array(df_all.loc[vinylpyrrolidoneb,'bio_rank_num3_b'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='D')
# polystyrene
Y = df_all.loc[polystyreneb,'LogP/SA (Å-2)'].values
X = df_all.loc[polystyreneb,'Tg (°C)'].values
plt_colors = df_all.loc[polystyreneb,'color_cat3_b'].tolist()
plt_areas = (np.array(df_all.loc[polystyreneb,'bio_rank_num3_b'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='P')
# pvc
Y = df_all.loc[pvcb,'LogP/SA (Å-2)'].values
X = df_all.loc[pvcb,'Tg (°C)'].values
plt_colors = df_all.loc[pvcb,'color_cat3_b'].tolist()
plt_areas = (np.array(df_all.loc[pvcb,'bio_rank_num3_b'].values) + 1)**1.5 * 25 
plt.scatter(X,Y,c=plt_colors,s=plt_areas,alpha=.5,linewidths=1,edgecolor='black',marker='d')

plt.ylabel('LogP/SA',fontsize=9)
plt.xlabel('Tg (C)',fontsize=9)
plt.tick_params(labelsize=9)
plt.ylim(-.01,.025)
plt.xlim(-75,200)
plt.tight_layout()
#plt.savefig(fname='Fig3a_062419.png',dpi=300,transparent=True,pad_inches=0.0)
plt.show()


# In[ ]:




