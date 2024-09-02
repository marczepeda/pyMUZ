### plot.py ###
# Author: Marc Zepeda
# Date: 2024-08-05

# Import packages
import os
import math
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Supporting methods
''' re_un_cap: Replace underscores with spaces and capitalizes each word for a given string
        input_string: input string
'''
def re_un_cap(input_string):
    output_string = input_string.replace('_', ' ').title()
    return output_string

''' round_up_pow_10: Rounds a given number to the nearest power of 10
        number: input number
'''
def round_up_pow_10(number):
    if number == 0:
        return 0

    exponent = math.ceil(math.log10(abs(number)))
    rounded = math.ceil(number / 10 ** exponent) * 10 ** exponent
    return rounded

''' log10: Returns log10 of maximum value from series or 0
        series: series with values
'''
def log10(series):
    return np.log10(np.maximum(series, 1))

''' move_dis_legend: Moves legend for distribution graphs.
        ax: matplotlib axis
        legend loc: legend location
'''
def move_dist_legend(ax,legend_loc,legend_title_size,legend_size,legend_bbox_to_anchor,legend_ncol):
    old_legend = ax.legend_
    handles = old_legend.legendHandles
    labels = [t.get_text() for t in old_legend.get_texts()]
    title = old_legend.get_title().get_text()
    ax.legend(handles,labels,loc=legend_loc,bbox_to_anchor=legend_bbox_to_anchor,
              title=title,title_fontsize=legend_title_size,fontsize=legend_size,ncol=legend_ncol)

''' extract_pivots: Returns a dictionary of pivot-formatted dataframes from tidy-formatted dataframe
        df: Tidy-formatted dataframe
        inds: y-axis
        cols: x-axis
        vars: variable column name (optional)
        vals: value column name (optional)
    Dependencies: pandas
'''
def extract_pivots(df: pd.DataFrame, x: str, y: str, vars='variable', vals='value'):
    piv_keys = list(df[vars].value_counts().keys())
    pivots = dict()
    for key in piv_keys:
        pivots[key]=pd.pivot(df[df[vars]==key],index=y,columns=x,values=vals)
    return pivots

''' formatter: Formats, displays, and saves plots.
        typ: plot type
        ax: matplotlib subplot
        df: tidy dataframe
    Dependencies: os, matplotlib, seaborn, re_un_cap(), round_up_pow_10():
'''
def formatter(typ,ax,df,x='x',y='y',cols='cols',file=None,dir=None,color_palette='colorblind',
              title='',title_size=18,title_weight='bold',
              x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_scale='linear',x_axis_dims=(0,100),x_ticks_rot=0,xticks=[],
              y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_scale='linear',y_axis_dims=(0,100),y_ticks_rot=0,yticks=[],
              legend_title='',legend_title_size=12,legend_size=12,legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_items=(0,0),legend_ncol=1):
    
    # Define plot types
    scats = ['scat', 'line', 'line_scat']
    cats = ['bar', 'box', 'violin', 'swarm', 'strip', 'point', 'count', 'bar_strip', 'box_strip', 'violin_strip','bar_swarm', 'box_swarm', 'violin_swarm']
    dists = ['hist', 'kde', 'hist_kde','rid']
    heats = ['ht']

    # Set color scheme
    sns.color_palette(color_palette)

    if typ not in heats:
        # Set title
        if title=='' and file is not None: title=re_un_cap(file[-4])
        plt.title(title, fontsize=title_size, fontweight=title_weight)
        
        # Set x axis
        if x_axis=='': x_axis=re_un_cap(x)
        plt.xlabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight)
        if typ not in cats: 
            plt.xscale(x_axis_scale)
            if x_axis_dims==(0,100) and x_axis_scale=='log': plt.xlim(0,round_up_pow_10(max(df[x])))
            else: plt.xlim(x_axis_dims[0],x_axis_dims[1])
        if xticks==[]: 
            if (x_ticks_rot==0)|(x_ticks_rot==90): plt.xticks(rotation=x_ticks_rot,ha='center')
            else: plt.xticks(rotation=x_ticks_rot,ha='right')
        else: 
            if (x_ticks_rot==0)|(x_ticks_rot==90): plt.xticks(ticks=xticks,labels=xticks,rotation=x_ticks_rot, ha='center')
            else: plt.xticks(ticks=xticks,labels=xticks,rotation=x_ticks_rot,ha='right')

        # Set y axis
        if y_axis=='': y_axis=re_un_cap(y)
        plt.ylabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight)
        plt.yscale(y_axis_scale)
        if y_axis_dims==(0,100) and y_axis_scale=='log': plt.ylim(0,round_up_pow_10(max(df[y])))
        else: plt.ylim(y_axis_dims[0],y_axis_dims[1])
        if yticks==[]: plt.yticks(rotation=y_ticks_rot)
        else: plt.yticks(ticks=yticks,labels=yticks,rotation=y_ticks_rot)

        # Set legend
        if legend_title=='': legend_title=cols
        if legend_items==(0,0) and typ not in dists:
            ax.legend(title=legend_title,title_fontsize=legend_title_size,fontsize=legend_size,
                    bbox_to_anchor=legend_bbox_to_anchor,loc=legend_loc,ncol=legend_ncol) # Move legend to the right of the graph
        elif typ not in dists:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(title=legend_title,title_fontsize=legend_title_size,fontsize=legend_size,
                    bbox_to_anchor=legend_bbox_to_anchor,loc=legend_loc,ncol=legend_ncol, # Move right of the graph
                    handles=handles[legend_items[0]:legend_items[1]],labels=labels[legend_items[0]:legend_items[1]]) # Only retains specified labels
        else: move_dist_legend(ax,legend_loc,legend_title_size,legend_size,legend_bbox_to_anchor,legend_ncol)

    # Save & show fig
    if file is not None and dir is not None:
        if not os.path.exists(dir):
            os.mkdir(dir)
        plt.savefig(fname=os.path.join(dir, file), dpi=600, bbox_inches='tight')
    plt.show()

# Graph methods
''' scat: Creates scatter plot related graphs.
        typ: plot type (scat, line, line_scat)
        df: tidy dataframe
        x: x-axis column
        y: y-axis column
    Dependencies: os, matplotlib, seaborn, formatter(), re_un_cap(), round_up_pow_10()
'''
def scat(typ: str,df: pd.DataFrame,x: str,y: str,cols=None,cols_ord=None,stys=None,cols_exclude=None,
         file=None,dir=None,color_palette='colorblind',edgecol='black',
         figsize=(10,6),title='',title_size=18,title_weight='bold',
         x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_scale='linear',x_axis_dims=(0,100),x_ticks_rot=0,xticks=[],
         y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_scale='linear',y_axis_dims=(0,100),y_ticks_rot=0,yticks=[],
         legend_title='',legend_title_size=12,legend_size=9,legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_items=(0,0),legend_ncol=1,
         **kwargs):

    # Omit excluded data
    if type(cols_exclude)==list: 
        for exclude in cols_exclude: df=df[df[cols]!=exclude]
    else: df=df[df[cols]!=cols_exclude]

    fig, ax = plt.subplots(figsize=figsize)
    
    if typ=='scat': sns.scatterplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, style=stys, edgecolor=edgecol, ax=ax, **kwargs)
    elif typ=='line': sns.lineplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, style=stys, ax=ax, **kwargs)
    elif typ=='line_scat':
        sns.lineplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, style=stys, ax=ax, **kwargs)  
        sns.scatterplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, style=stys, edgecolor=edgecol, ax=ax, **kwargs)
    else:
        print("Invalid type! scat, line, or line_scat")
        return
    
    formatter(typ,ax,df,x,y,cols,file,dir,color_palette,
              title,title_size,title_weight,
              x_axis,x_axis_size,x_axis_weight,x_axis_scale,x_axis_dims,x_ticks_rot,xticks,
              y_axis,y_axis_size,y_axis_weight,y_axis_scale,y_axis_dims,y_ticks_rot,yticks,
              legend_title,legend_title_size,legend_size,legend_bbox_to_anchor,legend_loc,legend_items,legend_ncol)

''' cat: Creates category dependent graphs.
        typ: plot type (bar, box, violin, swarm, strip, point, count, bar_swarm, box_swarm, violin_swarm)
        df: tidy dataframe
        x: x-axis column
        y: y-axis column
    Dependencies: os, matplotlib, seaborn, formatter(), re_un_cap(), round_up_pow_10()
'''
def cat(typ: str,df: pd.DataFrame,x: str,y: str,errorbar=None,cols=None,cols_ord=None,cols_exclude=None,
        file=None,dir=None,color_palette='colorblind',edgecol='black',lw=1,
        figsize=(10,6),title='',title_size=18,title_weight='bold',
        x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_scale='linear',x_axis_dims=(0,100),x_ticks_rot=0,xticks=[],
        y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_scale='linear',y_axis_dims=(0,100),y_ticks_rot=0,yticks=[],
        legend_title='',legend_title_size=12,legend_size=9,legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_items=(0,0),legend_ncol=1, 
        **kwargs):
    
    # Omit excluded data
    if type(cols_exclude)==list: 
        for exclude in cols_exclude: df=df[df[cols]!=exclude]
    else: df=df[df[cols]!=cols_exclude]

    fig, ax = plt.subplots(figsize=figsize)

    if typ=='bar': sns.barplot(data=df, x=x, y=y, errorbar=errorbar, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
    elif typ=='box': sns.boxplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, linewidth=lw, ax=ax, **kwargs)
    elif typ=='violin': sns.violinplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord,  edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
    elif typ=='swarm': sns.swarmplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, ax=ax, **kwargs)
    elif typ=='strip': 
        if cols is None: sns.stripplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
        else: sns.stripplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, ax=ax, **kwargs)
    elif typ=='point': sns.pointplot(data=df, x=x, y=y, errorbar=errorbar, hue=cols, hue_order=cols_ord, ax=ax, **kwargs)
    elif typ=='count': sns.countplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, ax=ax, **kwargs)
    elif typ=='bar_strip':
        sns.barplot(data=df, x=x, y=y, errorbar=errorbar, hue=cols, hue_order=cols_ord, edgecolor=edgecol, errcolor=edgecol, linewidth=lw, ax=ax, **kwargs)
        if cols is None:  sns.stripplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
        else: sns.stripplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, ax=ax, **kwargs)
    elif typ=='box_strip':
        sns.boxplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, linewidth=lw, ax=ax, **kwargs)
        if cols is None:  sns.stripplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
        else: sns.stripplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, ax=ax, **kwargs)
    elif typ=='violin_strip':
        sns.violinplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
        if cols is None:  sns.stripplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
        else: sns.stripplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, ax=ax, **kwargs)
    elif typ=='bar_swarm':
        sns.barplot(data=df, x=x, y=y, errorbar=errorbar, hue=cols, hue_order=cols_ord, edgecolor=edgecol, errcolor=edgecol, linewidth=lw, ax=ax, **kwargs)
        if cols is None:  sns.swarmplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
        else: sns.swarmplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, ax=ax, **kwargs)
    elif typ=='box_swarm':
        sns.boxplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, linewidth=lw, ax=ax, **kwargs)
        if cols is None:  sns.swarmplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
        else: sns.swarmplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, ax=ax, **kwargs)
    elif typ=='violin_swarm':
        sns.violinplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
        if cols is None:  sns.swarmplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
        else: sns.swarmplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, ax=ax, **kwargs)
    else:
        print('Invalid type! bar, box, violin, swarm, strip, point, count, bar_strip, box_strip, violin_strip, bar_swarm, box_swarm, violin_swarm')
        return

    formatter(typ,ax,df,x,y,cols,file,dir,color_palette,
              title,title_size,title_weight,
              x_axis,x_axis_size,x_axis_weight,x_axis_scale,x_axis_dims,x_ticks_rot,xticks,
              y_axis,y_axis_size,y_axis_weight,y_axis_scale,y_axis_dims,y_ticks_rot,yticks,
              legend_title,legend_title_size,legend_size,legend_bbox_to_anchor,legend_loc,legend_items,legend_ncol)

''' dist: Creates distribution graphs.
        typ: plot type (hist, kde, hist_kde, rid)
        df: tidy dataframe
        x: x-axis column
    Dependencies: os, matplotlib, seaborn, formatter(), re_un_cap(), round_up_pow_10()
'''
def dist(typ: str,df: pd.DataFrame,x: str,cols=None,cols_ord=None,cols_exclude=None,bins=40,log10_low=0,
        file=None,dir=None,color_palette='colorblind',edgecol='black',lw=1,ht=1.5,asp=5,tp=.8,hs=0,des=False,
        figsize=(10,6),title='',title_size=18,title_weight='bold',
        x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_scale='linear',x_axis_dims=(0,100),x_ticks_rot=0,xticks=[],
        y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_scale='linear',y_axis_dims=(0,100),y_ticks_rot=0,yticks=[],
        legend_title='',legend_title_size=12,legend_size=9,legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_items=(0,0),legend_ncol=1, 
        **kwargs):

    # Omit excluded data
    if type(cols_exclude)==list: 
        for exclude in cols_exclude: df=df[df[cols]!=exclude]
    else: df=df[df[cols]!=cols_exclude]

    if typ=='hist':
        fig, ax = plt.subplots(figsize=figsize)
        if isinstance(bins, int):
            if x_axis_scale=='log':
                bins = np.logspace(log10(df[x]).min(), log10(df[x]).max(), bins + 1)
            else:
                bins = np.linspace(df[x].min(), df[x].max(), bins + 1)
        sns.histplot(data=df, x=x, kde=False, bins=bins, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
        y='Count'
        formatter(typ,ax,df,x,y,cols,file,dir,color_palette,
                  title,title_size,title_weight,
                  x_axis,x_axis_size,x_axis_weight,x_axis_scale,x_axis_dims,x_ticks_rot,xticks,
                  y_axis,y_axis_size,y_axis_weight,y_axis_scale,y_axis_dims,y_ticks_rot,yticks,
                  legend_title,legend_title_size,legend_size,legend_bbox_to_anchor,legend_loc,legend_items,legend_ncol)
    elif typ=='kde': 
        fig, ax = plt.subplots(figsize=figsize)
        if x_axis_scale=='log':
            df[f'log10({x})']=np.maximum(np.log10(df[x]),log10_low)
            df=df.reset_index(drop=True)
            sns.kdeplot(data=df, x=f'log10({x})', hue=cols, hue_order=cols_ord, linewidth=lw, ax=ax, **kwargs)
            x_axis_scale='linear'
            if x_axis=='': x_axis=f'log10({x})'
        else: sns.kdeplot(data=df, x=x, hue=cols, hue_order=cols_ord, linewidth=lw, ax=ax, **kwargs)
        y='Density'
        formatter(typ,ax,df,x,y,cols,file,dir,color_palette,
                  title,title_size,title_weight,
                  x_axis,x_axis_size,x_axis_weight,x_axis_scale,x_axis_dims,x_ticks_rot,xticks,
                  y_axis,y_axis_size,y_axis_weight,y_axis_scale,y_axis_dims,y_ticks_rot,yticks,
                  legend_title,legend_title_size,legend_size,legend_bbox_to_anchor,legend_loc,legend_items,legend_ncol)
    elif typ=='hist_kde':
        fig, ax = plt.subplots(figsize=figsize)
        if x_axis_scale=='log':
            df[f'log10({x})']=np.maximum(np.log10(df[x]),log10_low)
            df=df.reset_index(drop=True)
            sns.histplot(data=df, x=f'log10({x})', kde=True, bins=bins, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
            x_axis_scale='linear'
            if x_axis=='': x_axis=f'log10({x})'
        else: sns.histplot(data=df, x=x, kde=True, bins=bins, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
        y='Count'
        formatter(typ,ax,df,x,y,cols,file,dir,color_palette,
                  title,title_size,title_weight,
                  x_axis,x_axis_size,x_axis_weight,x_axis_scale,x_axis_dims,x_ticks_rot,xticks,
                  y_axis,y_axis_size,y_axis_weight,y_axis_scale,y_axis_dims,y_ticks_rot,yticks,
                  legend_title,legend_title_size,legend_size,legend_bbox_to_anchor,legend_loc,legend_items,legend_ncol)
    elif typ=='rid':
        sns.color_palette(color_palette)
        if x_axis_scale=='log':
            df[f'log10({x})']=np.maximum(np.log10(df[x]),log10_low)
            df=df.reset_index(drop=True)
            g = sns.FacetGrid(df, row=cols, hue=cols, col_order=cols_ord, hue_order=cols_ord, height=ht, aspect=asp)
            g.map(sns.kdeplot, f'log10({x})', linewidth=lw, **kwargs)
            if x_axis=='': x_axis=f'log10({x})'
        else:
            g = sns.FacetGrid(df, row=cols, hue=cols, col_order=cols_ord, hue_order=cols_ord, height=ht, aspect=asp)
            g.map(sns.kdeplot, x, linewidth=lw, **kwargs)
            if x_axis=='': x_axis=x
        for ax in g.axes.flatten():
            ax.set_xlim(x_axis_dims[0],x_axis_dims[1])
            ax.set_xlabel(x_axis,fontsize=x_axis_size,fontweight=x_axis_weight)
        g.set(yticks=yticks, ylabel=y_axis)
        g.set_titles("")
        if title=='' and file is not None: title=re_un_cap(file[-4])
        g.figure.suptitle(title, fontsize=title_size, fontweight=title_weight)
        g.figure.subplots_adjust(top=tp,hspace=hs)
        if des==False: g.despine(top=False,right=False)
        else: g.despine(left=True)
        if legend_title=='': legend_title=cols
        g.figure.legend(title=legend_title,title_fontsize=legend_title_size,fontsize=legend_size,
                        loc=legend_loc,bbox_to_anchor=legend_bbox_to_anchor)
        if file is not None and dir is not None:
            if not os.path.exists(dir):
                os.mkdir(dir)
            plt.savefig(fname=os.path.join(dir, file), dpi=600, bbox_inches='tight')
        plt.show()
    else:
        print('Invalid type! hist, kde, hist_kde, rid')
        return
        
''' heat: Creates heat plot related graphs.
        df: tidy dataframe
        x: x-axis
        y: y-axis
    Dependencies: os, matplotlib, seaborn, formatter(), re_un_cap(), round_up_pow_10()
'''
def heat(df: pd.DataFrame, x: str, y: str, vars='variable', vals='value', typ='ht',
         file=None,dir=None,color_palette='colorblind',edgecol='black',lw=1,annot=False,cmap="Reds",sq=True,
         title='',title_size=18,title_weight='bold',
         x_axis='',x_axis_size=12,x_axis_weight='bold',x_ticks_rot=0,
         y_axis='',y_axis_size=12,y_axis_weight='bold',y_ticks_rot=0,
         **kwargs):

    # Find min and max values in the dataset for normalization
    vmin = df[vals].values.min()
    vmax = df[vals].values.max()

    # Generates dictionary of pivot tables
    pivots = extract_pivots(df=df,x=x,y=y,vars=vars,vals=vals)
    pivots['colorbar']=pivots[next(iter(pivots))]

    # Create a single figure with multiple heatmap subplots
    fig, axes = plt.subplots(nrows=1,ncols=len(list(pivots.keys())),figsize=(pivots[next(iter(pivots))].shape[0]*len(list(pivots.keys())),pivots[next(iter(pivots))].shape[1]),sharex=True,sharey=True)
    for i, (ax, subtitle) in enumerate(zip(axes, list(pivots.keys()))):
        if i!=len(axes)-1: sns.heatmap(pivots[subtitle],annot=annot,cmap=cmap,square=sq,ax=ax,vmin=vmin,vmax=vmax,linecolor=edgecol,linewidths=lw,cbar=False) # No heatmap colorbar
        else: # Add heatmap colorbar to the end
            sns.heatmap(pivots[subtitle],annot=annot,cmap=cmap,ax=ax,vmin=vmin,vmax=vmax,linecolor=edgecol,linewidths=lw,cbar=True)
            plt.gca().set_visible(False)
        if len(list(pivots.keys()))>1: ax.set_title(subtitle,fontsize=x_axis_size,fontweight=x_axis_weight)  # Add title to subplot
        else: ax.set_title(title,fontsize=x_axis_size,fontweight=x_axis_weight)
        if i==0: # First heatmap subplot has x & y axes
            if x_axis=='': ax.set_xlabel(re_un_cap(x),fontsize=x_axis_size,fontweight=x_axis_weight)
            else: ax.set_xlabel(x_axis,fontsize=x_axis_size,fontweight=x_axis_weight)
            if y_axis=='': ax.set_ylabel(re_un_cap(y),fontsize=y_axis_size,fontweight=y_axis_weight)
            else: ax.set_ylabel(y_axis,fontsize=y_axis_size,fontweight=y_axis_weight)
            ax.tick_params(axis='x',rotation=x_ticks_rot)
            ax.tick_params(axis='y',rotation=y_ticks_rot)
        else: # Remove redundant information from remaining heatmap subplots 
            if x_axis=='': ax.set_xlabel(re_un_cap(x),fontsize=x_axis_size,fontweight=x_axis_weight)
            else: ax.set_xlabel(x_axis,fontsize=x_axis_size,fontweight=x_axis_weight)
            ax.set_ylabel('')
            ax.tick_params(axis='x',rotation=x_ticks_rot)
        ax.set_facecolor('white')  # Set background to transparent

    formatter(typ,ax,df,x,y,'cols',file,dir,color_palette,
              title,title_size,title_weight,
              x_axis,x_axis_size,x_axis_weight,'linear',(0,100),x_ticks_rot,
              y_axis,y_axis_size,y_axis_weight,'linear',(0,100),y_ticks_rot,
              '',12,9,(1,1),'upper left',(0,0))
    
''' stack: Creates stacked bar plot
        df: dataframe from outcomes
        x: x-axis column
        y: y-axis column
        cols: colors column
        cutoff: y-axis values needs be greater than (ex: 0%)
    Dependencies: plot.py,re,os,pandas,numpy,matplotlib.pyplot
'''
def stack(df: pd.DataFrame,x='sample',y='fraction',cols='edit',cutoff=0,cols_ord=[],
          file=None,dir=None,color_palette='viridis',
          figsize=(10,6),title='',title_size=18,title_weight='bold',
          x_axis='',x_axis_size=12,x_axis_weight='bold',x_ticks_rot=45,
          y_axis='',y_axis_size=12,y_axis_weight='bold',y_ticks_rot=0,
          legend_title='',legend_title_size=12,legend_size=12,
          legend_bbox_to_anchor=(1,1),legend_loc='upper left',**kwargs):
    
    # Make pivot table
    df_cut=df[df[y]>cutoff]
    df_pivot=pd.pivot_table(df_cut, index=x, columns=cols, values=y, aggfunc='mean')
    if cols_ord!=[]: df_pivot=df_pivot[cols_ord]
    
    # Make stacked barplot
    df_pivot.plot(kind='bar',figsize=figsize,colormap=color_palette,stacked=True,**kwargs)

    # Set title
    if title=='' and file is not None: title=re_un_cap(file[-4])
    plt.title(title, fontsize=title_size, fontweight=title_weight)
    
    # Set x axis
    if x_axis=='': x_axis=re_un_cap(x)
    plt.xlabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight)
    if x_ticks_rot!=90: plt.xticks(rotation=x_ticks_rot, ha='right')

    # Set y axis
    if y_axis=='': y_axis=re_un_cap(y)
    plt.ylabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight)
    plt.yticks(rotation=y_ticks_rot)

    # Set legend
    if legend_title=='': legend_title=re_un_cap(cols)
    plt.legend(title=legend_title, title_fontsize=legend_title_size, fontsize=legend_size, 
               bbox_to_anchor=legend_bbox_to_anchor, loc=legend_loc)
    
    # Save & show fig
    if file is not None and dir is not None:
        if not os.path.exists(dir):
            os.mkdir(dir)
        plt.savefig(fname=os.path.join(dir, file), dpi=600, bbox_inches='tight')
    plt.show()