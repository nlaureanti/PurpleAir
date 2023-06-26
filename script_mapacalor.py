#!/home/nicole/miniconda3/envs/shapefile_cut/bin/python
# coding: utf-8

# In[1]:
# conda activate Purpleair (?)


import pandas as pd    # conjunto de funções que se usa para dados em tabela (csv)
#https://pandas.pydata.org/docs/index.html
import matplotlib.pyplot as plt # conjunto de funções que se usa para plotagem 
#https://matplotlib.org/stable/index.html
import numpy as np #para trabalhar com vetores/arrays


import cartopy.crs as ccrs #para introduzir os eixos cartográficos

from scipy.interpolate import griddata as gd #para interpolar os dados
from mpl_toolkits.axes_grid1.inset_locator import inset_axes #para movimentar o eixo que possui a colorbar
from scipy.ndimage.filters import gaussian_filter  #filtro para suavizar o campo pós-interpolação
import matplotlib.transforms as mtransforms # para transformar eixos


# # Fazendo a leitura das coordenadas

# In[2]:


fname_locais='Números_ID_PA.xlsx'
df_locais=pd.read_excel(fname_locais,
                      header=2, index_col=0, usecols=["LOCAL","COORD LAT","COORD LONG"])
print(df_locais)


# # Arquivos de entrada (preparando a execução):

# In[3]:


fnames=['Ilha do Mel corrigido.csv',        'Reserva dos Papagaios corrigido.csv',
'Reserva das Águas corrigido.csv',              'RF 20211010-20211016 corrigido.csv',
'Reserva da Guaricica (outside) (-25.314958 -48.696156) Primary Real Time 05_05_2022 05_07_2022-corrigido.csv'
]
pasta_arq='./corrigidos_mapadecalor/'


# In[4]:


#Atribuir a coordenada ao nome do arquivo
vetor_lat=[] ;vetor_lon=[]
list_estacao=[]
name_estacao=[]
for loc in range(len(df_locais.index)):
    for f in fnames:
        if df_locais.index[loc] in f.upper(): 
            #condição para verificação do nome do arquivo e atribuição de coordenada
            print(f"loc: {loc} {df_locais.index[loc]} ok {f.upper()}")
            vetor_lat.append(df_locais["COORD LAT"][loc])
            vetor_lon.append(df_locais["COORD LONG"][loc])
            list_estacao.append(f"{pasta_arq}/{f}")
            name_estacao.append(df_locais.index[loc])
print(name_estacao) 
print(vetor_lat)
print(vetor_lon) 


# In[5]:


#função para fazer a abertura dos arquivos
def get_data(f):
    try:
        df=pd.read_csv(f,
                      sep=',',header=0, index_col=0,
                       parse_dates=["dates"], na_values='NaN') #arquivo online
    except:
        df=pd.read_csv(f,
                      sep=',',header=0, index_col=0,
                       parse_dates=["UTCDateTime"], na_values='NaN') #arquivo do SD
    return df


# In[6]:


#Crio o grid 
dfs=[];xs=[];ys=[]
for i in range(len(name_estacao)):
    print(name_estacao[i])
    xs.append(vetor_lon[i])
    ys.append(vetor_lat[i])
    
    dfs.append(get_data(list_estacao[i]))
    print(xs[i],ys[i])
    #print(loc,xs,ys,name_estacao[i],'\n',dfs[i])
xs, ys = np.array(xs), np.array(ys)  #criar o grid dos dados


# In[7]:


#função para interpolar o dado dentro do grid
def interpolar_xy(datax,datay,dataz):
    # Interpolation: Generate grid data of latitude and longitude
    numcols, numrows = 100, 100
    xi = np.linspace(np.min(datax), np.max(datax), numcols)
    yi = np.linspace(np.min(datay), np.max(datay), numrows)
    xi, yi = np.meshgrid(xi, yi)    
     
    zi = gd(
            (datax,datay),
            dataz,
            (xi, yi),
            method='nearest')
    zi = gaussian_filter(zi, 5) 
    return xi, yi, zi    


# In[8]:


#Função para plotar o arquivo shapefile
def plot_shape(fname,ax,color='none'):
    from cartopy.io.shapereader import Reader
    """ PLOTS EVERY SHAPE """

    ax.add_geometries(Reader(fname).geometries(),
                  ccrs.PlateCarree(),
                  facecolor=color, hatch=None, edgecolor='black')


# # 1º Plot: mapa da estações

# In[9]:


#[ print(f" {name_estacao[i]} {xs[i]} {ys[i]} \n {matrix[i]} ") for i in range(len(name_estacao)) ]

proj = ccrs.PlateCarree() 
#cria a projeção do mapa # https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html
fig, ax = plt.subplots(1,figsize=(15,15),
                   subplot_kw=dict(projection=proj),
                   sharex=True,sharey=True,
                   gridspec_kw={'hspace': 0, 'wspace':0}
                      )
for d in range(1):
    ax.set_extent((-49.0, -48.0, -25.80, -25.), crs=proj) #define os limites do mapa
    shpname='shoreNOAA_GSHHG_cutPR_new.shp'
    plot_shape(shpname,ax,'goldenrod')   #chamando a função para plotar o shapefile
    
    [ ax.text(xs[i],ys[i],f"{name_estacao[i]}",
              color='black',fontsize=16) for i in range(len(name_estacao)) ] #plotar os nomes das estações
    [ ax.plot(xs[i],ys[i],color='white',marker='^',
                    markerfacecolor='red',
              markersize=9, transform=proj) for i in range(len(name_estacao)) ] #plotar os marcadores
    #https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html

    #Linhas de grade para definir o domínio em latitude e longitude
    gl = ax.gridlines(crs=proj, draw_labels=True,
                                linewidth=1, color='k',
                                alpha=0.2, linestyle='--')        
    gl.right_labels = None
    gl.top_labels = None
    
    #Criar o box com o nome "Estações"
    trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)            
    ax.text(0.0, 1.0, "Estações", transform=ax.transAxes + trans,
            fontsize=20, verticalalignment='top', 
            bbox=dict(facecolor='1', edgecolor='k', pad=3.0))
    
    plt.savefig('mapa_Estacoes.png',facecolor='white')
    print(' > mapa_Estacoes.png')


# # Leitura e organização dos dados 

# In[10]:


df_horario=[]
matrix_hourlymean=[]; matrix_hourlymin=[]; matrix_hourlymax=[]
matrix_dailymean=[]; matrix_dailymin=[]; matrix_dailymax=[]
for i in range(len(name_estacao)):
    datai=dfs[i].index[0] #utilize esses dois parmâetros para alterar a data de leitura
    dataf=dfs[i].index[-1]
    date_range=pd.date_range(start=datai,
                    end=dataf,
                    freq="h") #controla o nome dos arquivos/planilhas

    print(f"usando resample em {name_estacao[i]} para {dfs[i].resample('60T')}")

    matrix_hourlymean.append(dfs[i].groupby(dfs[i].index.hour)['pm2.5'].mean())
    matrix_hourlymax.append(dfs[i].groupby(dfs[i].index.hour)['pm2.5'].max())
    matrix_hourlymin.append(dfs[i].groupby(dfs[i].index.hour)['pm2.5'].min())
    
    matrix_dailymean.append(dfs[i].resample('D')['pm2.5'].mean())
    matrix_dailymax.append(dfs[i].resample('D')['pm2.5'].max())
    matrix_dailymin.append(dfs[i].resample('D')['pm2.5'].min())
    #https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.min.html


# In[11]:


#Saida para verificação com as planilhas originais
for matrix in [matrix_hourlymean, matrix_dailymean]:
    [ print(f" {name_estacao[i]} {xs[i]} {ys[i]} \n {matrix[i]} ") for i in range(len(name_estacao)) ]


# In[24]:


#Função para plotar o mapa de calor
def plotar_mapa_calor(matriz, nt, c_levels=np.linspace(0,20,11), titulo_str= 'PM2.5', dayformat='%d/%m/%y'):
    proj = ccrs.PlateCarree() #projeção
    if nt<10:
        #quando o número de gáficos é menor do que 10: plota em uma única coluna
        fig, ax = plt.subplots(nt,figsize=(29,34),
                       subplot_kw=dict(projection=proj),
                       sharex=True,sharey=True,
                       gridspec_kw={'hspace': 0, 'wspace':0}
                          )
        for d in range(nt):
            drange=matriz[0].index[d]
            #ax[d].set_extent((-49.0, -48.0, -25.80, -25.), crs=proj)
            shpname='shoreNOAA_GSHHG_cutPR.shp'
            plot_shape(shpname,ax[d])  #plotar shapefile

            az=np.empty(len(name_estacao)); #criar matriz para interpolação
            for i in range(len(name_estacao)):
                try:
                    az[i] = matriz[i][d] #recebe dado de cada estação i para o dia d
                except:
                    az[i]=np.nan
    #        az=np.array([matriz[i][d] for i in range(len(name_estacao))])
            print(f"day {d} {drange}",xs.shape,ys.shape,az.shape)
            x,y,z=interpolar_xy(xs,ys,az) #interpola
            
          
            #plotar os valores no gráfico?
            [ ax[d].plot(xs[i],ys[i],label=name_estacao,color='white',marker='o',
                            markerfacecolor='black',markersize=4, transform=proj)
                         for i in range(len(name_estacao)) ] 
            [ ax[d].text(xs[i],ys[i],f"{az[i]:.2f}") for i in range(len(name_estacao)) ]
            
            #Plotar o mapa interpolado:
            p=ax[d].contourf(x,y,z, levels=c_levels, cmap='hot_r', extend='max')    


            
            #ajuste dos grids de coordenadas 
            gl = p.axes.gridlines(crs=proj, draw_labels=True,
                                        linewidth=1, color='k',
                                        alpha=0.2, linestyle='--')        
            gl.right_labels = None
            gl.top_labels = None     
            if d != nt-1:
                gl.bottom_labels = None

            #título plotado no canto superior à esquerda: a data, nesse caso
            trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)            
            ax[d].text(0.0, 1.0, pd.to_datetime(drange).strftime(dayformat), transform=ax[d].transAxes + trans,
                    fontsize='medium', verticalalignment='top', 
                    bbox=dict(facecolor='1', edgecolor='k', pad=3.0))

        #ax[0].set_title(pd.to_datetime(drange).strftime('%D-%M-%Y'),loc='left')
        #colobar:
        cbar_ax = inset_axes(ax[d],
                        width="3%",  
                        height="72%",
                        loc='center right',
                        borderpad=-4,
                       )
        cb = fig.colorbar(p, orientation='vertical', label='$mu$/hora',
                            ticks=c_levels,
                           aspect=30,shrink=.6,pad=3,cax=cbar_ax)
        cb.ax.tick_params(labelsize=15)
        
        #titulo centralizado superior:
        fig.suptitle(titulo_str
                ,fontsize=20,y=.89)
        plt.savefig(f"{titulo_str}.png",facecolor='white')
        print(f" >{titulo_str}.png")
    else:
        fig, ax = plt.subplots(12,2,figsize=(29,34),
                       subplot_kw=dict(projection=proj),
                       sharex=True,sharey=True,
                       gridspec_kw={'hspace': 0, 'wspace':0}
                          )
        n=0;m=0
        for d in range(nt):
            drange=matriz[0].index[d]
            #ax[n,m].set_extent((-49.0, -48.0, -25.80, -25.), crs=proj) #para alterar os limites do mapa
            shpname='shoreNOAA_GSHHG_cutPR_new.shp'
            plot_shape(shpname,ax[n,m])               

            az=np.empty(len(name_estacao));
            for i in range(len(name_estacao)):
                try:
                    az[i] = matriz[i][d]
                except:
                    az[i]=np.nan
    #        az=np.array([matriz[i][d] for i in range(len(name_estacao))])
            print(f"day {d} {drange}",xs.shape,ys.shape,az.shape)
            x,y,z=interpolar_xy(xs,ys,az)
            
         

            [ ax[n,m].plot(xs[i],ys[i],label=name_estacao,color='white',marker='o',
                            markerfacecolor='black',markersize=4, transform=proj)
                     for i in range(len(name_estacao)) ] 
            [ ax[n,m].text(xs[i],ys[i],f"{az[i]:.2f}") for i in range(len(name_estacao)) ]
            
            p=ax[n,m].contourf(x,y,z, levels=c_levels, cmap='hot_r', extend='max')   
            #opções de cores https://matplotlib.org/stable/tutorials/colors/colormaps.html


            gl = p.axes.gridlines(crs=proj, draw_labels=True,
                                        linewidth=1, color='k',
                                        alpha=0.2, linestyle='--')        
            gl.right_labels = None
            gl.top_labels = None

            #ax[n,m].coastlines()

            

            trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)            
            ax[n,m].text(0.0, 1.0, pd.to_datetime(matriz[0].index[d]).strftime(dayformat), transform=ax[n,m].transAxes + trans,
                    fontsize='medium', verticalalignment='top', 
                    bbox=dict(facecolor='1', edgecolor='k', pad=3.0))
            n+=1
            if n!=nt/2:
                gl.bottom_labels = None
            if n>=nt/2: m+=1; n=0
            if m>=nt/2: break

        #ax[0].set_title(pd.to_datetime(drange).strftime('%D-%M-%Y'),loc='left')
        cbar_ax = inset_axes(ax[11,1],
                        width="3%",  
                        height="72%",
                        loc='center right',
                        borderpad=-4,
                       )
        cb = fig.colorbar(p, orientation='vertical', label='$mu$/hora',
                            ticks=c_levels,
                           aspect=30,shrink=.6,pad=3,cax=cbar_ax)
        cb.ax.tick_params(labelsize=15)
        fig.suptitle(titulo_str
                ,fontsize=20,y=.89)
        plt.savefig(f"{titulo_str}.png",facecolor='white')
        print(f" > {titulo_str}.png")  


# In[13]:



print(f"Matriz média diária {len(matrix_dailymean[0].index)}")

plotar_mapa_calor(matrix_dailymean, len(matrix_dailymean[0].index)-1, dayformat='%d/%m/%y',
                  c_levels=np.linspace(0,20,11), 
                  titulo_str= 'Daily Mean PM2.5')


# In[25]:


print(f"Matriz média horária {len(matrix_hourlymean[0].index)}")
plotar_mapa_calor(matrix_hourlymean, 24, dayformat= '%Hh',
                  c_levels=np.linspace(0,15,11), 
                  titulo_str= 'Hourly Mean PM2.5')


# In[19]:


print(f"Matriz max horária {len(matrix_hourlymax[0].index)}")
plotar_mapa_calor(matrix_hourlymax, 24, dayformat= '%Hhoras',
                  c_levels=np.linspace(0,250,11), 
                  titulo_str= 'Hourly Max PM2.5')


# In[16]:


print(f"Matriz máximo horário {len(matrix_dailymax[0].index)}")
plotar_mapa_calor(matrix_dailymax, len(matrix_dailymax[0].index), 
                  c_levels=np.linspace(0,250,11), dayformat='%d/%m/%Y',
                  titulo_str= 'Daily maximum PM2.5')





