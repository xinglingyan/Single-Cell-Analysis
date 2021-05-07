
import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib.pyplot as plt
plt.switch_backend('agg')
import seaborn as sns
sns.set_style( "white" )
plt.rcParams[ "font.size" ] = 4.0
plt.rcParams[ "figure.dpi" ] = 100
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['axes.labelsize'] = 6 
plt.rcParams[ "figure.figsize" ] = ( 2*0.8,2.75*0.8 )
plt.rcParams[ "font.serif" ] = 'Arial'

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80)

raw = sc.read_10x_mtx(
    './',  # the directory with the `.mtx` file
    var_names='gene_symbols',                  # use gene symbols for the variable names (variables-axis index)
    cache=True)    

adata.var_names_make_unique()
#filter gene that expressed in less than 10 cells
sc.pp.filter_genes(adata, min_cells=10)
#filter cell that less than 1000 genes
tmp = sc.pp.filter_cells(adata, min_genes=1000,inplace =False)
adata2 = adata[tmp[0],:]
# normalize data
sc.pp.normalize_total(adata2, target_sum=1e4)
sc.pp.log1p(adata2)
hv = sc.pp.highly_variable_genes(adata2,n_top_genes=4000,inplace=False)

adata3 = adata2[:, hv.highly_variable]
sc.pp.pca(adata3)

sc.pp.neighbors(adata3)
sc.tl.umap(adata3)
sc.pl.umap(adata3)
sc.tl.leiden(adata3,resolution=0.5,directed=False)
sc.pl.umap(adata3,color=['leiden','olig2','isl1'])
sc.pl.umap(adata3,color=['olig1',"olig2",'isl1','isl2a','lhx1a','sox10','nkx6.1', 'irx3a', 'LHX3','mnx2a','mnx2b','mnx1', 'gfap','fabp7a','vsx2','vsx1','gata3','tal1','gata2a','tal2','en1b'],save='markers.pdf')
sc.tl.rank_genes_groups(adata3, 'leiden', method='t-test',use_raw =False,n_genes =100)

marker_genes = {
    'pMN':{'olig1', 'olig2'},'Motor neuron':{'isl1', 'lhx1'},
    'OPC':{'olig1','olig2','sox10'}
   #'OPC':{'sox10'}
}

marker_matches = sc.tl.marker_gene_overlap(adata3,marker_genes,method ='overlap_count',adj_pval_threshold=0.01)
marker_matches

cluster_names = []
for col in marker_matches.columns:
    if marker_matches[col].sum() == 0:
        cluster_names.append('Unknown_{}'.format(col))
    else:
        name = marker_matches[col].sort_values().index[-1]
        cluster_names.append(name+'_'+col)

cluster_names = ['Motor neuron_0', 'Motor neuron_1', 'pMN_2', 'Radial glia_3', 'pMN_4', 'Motor neuron_5', 'Unknown_6', 'pMN > Motor neuron_7', 'Interneuron_8', 'Unknown_9', 'Interneuron_10', 'OPC_11', 'Radial_12', 'Unknown_13']
adata3.rename_categories('leiden', cluster_names)
sc.settings.set_figure_params(dpi=200,dpi_save=1000,frameon=True)
#save umap
sc.pl.umap(adata3,color=['olig1',"olig2",'isl1','lhx1a','sox10'],save='markers.pdf')
sc.pl.umap(adata3,color=['nkx6.1', 'irx3a','LHX3','mnx2a','mnx2b','mnx1', 'isl2a','isl2b'],save='markers_2.pdf')

# trajectory
sc.tl.paga(adata3, groups='leiden',model='v1.2')
sc.settings.set_figure_params(dpi=100,dpi_save=1000,frameon=True)
sc.pl.paga(adata3,color=['leiden'],save='abstract_paga.pdf')
sc.tl.draw_graph(adata3,init_pos='paga' ,layout ='fa')
sc.settings.set_figure_params(dpi=300,dpi_save=1000,frameon=True)
sc.pl.draw_graph(adata3, color=['leiden','sox10'], legend_loc='on data', title='', legend_fontsize=6,size=4,frameon =True,save='.pdf')
#pseudotime
x = adata3.copy()
sc.tl.diffmap(x)
sc.pp.neighbors(x, n_neighbors=10)
sc.tl.draw_graph(x)
sc.settings.set_figure_params(dpi=300,dpi_save=1000,frameon=True)
sc.pl.draw_graph(x, color=['leiden','sox10'], legend_loc='on data', title='', legend_fontsize=6,size=4,frameon =True)
adata3.uns['iroot'] = np.flatnonzero(adata3.obs['leiden']  == 'pMN_4')[0]
sc.tl.dpt(adata3)#, n_dcs=5,min_group_size =0.05,n_branchings =1)
sc.pl.draw_graph(adata3, color=['leiden','dpt_pseudotime'], legend_loc='on data', title='', 
    legend_fontsize=6,size=4,vmax = 0.5,color_map ='Blues',save='.pdf')
#save data
adata3.write('adata3.h5ad')
##########################################################################################
##########################################################################################
adata3.uns['leiden_colors']=['#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b',
       '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a',
       '#ff9896', '#c5b0d5']

cluster_names = ['Motor neuron_0', 'Motor neuron_1', 'pMN_2', 'Radial glia_3', 'pMN_4', 'Motor neuron_5', 'Unknown_6', 'pMN > Motor neuron_7', 'Interneuron_8', 'Unknown_9', 'Interneuron_10', 'OPC_11', 'Unknown_12', 'Unknown_13']
adata3.rename_categories('leiden', cluster_names)
sc.tl.paga(adata3, groups='leiden',model='v1.2')
adata3.uns['leiden_colors']=['#1f77b4', '#ff7f0e', '#279e68', 'lightgrey', '#aa40fc', '#8c564b',
       'lightgrey', '#b5bd61', 'lightgrey', 'lightgrey', 'lightgrey', '#98df8a','lightgrey', 'lightgrey']

adata3.uns['iroot'] = np.flatnonzero(adata3.obs['leiden']  == 'pMN_4')[0]
sc.tl.dpt(adata3, n_dcs=5,min_group_size =0.05,n_branchings =2)

sc.tl.rank_genes_groups(adata3,groupby='leiden', groups = ['OPC_11'],reference ='pMN > Motor neuron_7', method='t-test',use_raw=False,n_genes =10000 )



sc.settings.set_figure_params(dpi=100,dpi_save=1000,frameon=True)
sc.tl.rank_genes_groups(st_adata2,groupby='leiden', groups = ['OPC_11'],reference ='pMN > Motor neuron_7', method='t-test',n_genes =20000 )

x = sc.get.rank_genes_groups_df(st_adata2, group="OPC_11"); x.to_csv('pMNMotor7_OPC11.txt',sep='\t')
def get_TF_anno():
    f = 'Danio_rerio_TF.txt'
    tfs = {}
    for line in open(f):
        line = line.strip().split('\t')
        tfs[line[1]] = 'TF'
    return tfs

tfs = get_TF_anno()

def get_DEG(f):
    #f = 'pMNMotor7_OPC11.txt'
    df = pd.read_table(f,index_col=0);df.index = df['names']
    #df = df.sort_values(by='pvals_adj')
    ups = df[(df['pvals_adj']<0.001)&(df['logfoldchanges']>1)].index.tolist()
    print(len(ups))
    downs = df[(df['pvals_adj']<0.001)&(df['logfoldchanges']<-1)].index.tolist()
    print(len(downs))
    np.random.shuffle(ups);np.random.shuffle(downs)
    degs = ups + downs
    fo = f.replace('.txt','_DEGs.txt')
    with open(fo,'w') as o:
        for gene in degs:
            o.write(gene+'\n')
    return degs

degs_7_11 = get_DEG('pMNMotor7_OPC11.txt')


cells = st_adata2.obs
group_7_11_cells = cells[(cells['leiden']=='pMN > Motor neuron_7')|(cells['leiden']=='OPC_11')].index
group_7_11_st_adata2 = st_adata2[group_7_11_cells,:]

sc.pp.scale(group_7_11_st_adata2)
sc.pl.heatmap(group_7_11_st_adata2, degs_7_11, groupby='leiden', figsize=(5,7),
              #var_group_positions=[(0,10), (11, 12)], 
              vmin=-4, 
              vmax=4,
              cmap='vlag',
              var_group_rotation=0, dendrogram=False,
             #standard_scale ='var',
              swap_axes =True,
              show_gene_labels = False,
              save = 'heatmap_7_11_scale.pdf'
             )

def get_tf(f,degs):
    fo = f.replace('.txt','_TFs.txt')
    degs_tf = []
    for gene in degs:
        if gene in tfs:
            degs_tf.append(gene)
    with open(fo,'w') as o:
        for gene in degs_tf:
            o.write(gene+'\n')
    return degs_tf


degs_7_11_tf = get_tf('pMNMotor7_OPC11.txt',degs_7_11)
sc.pl.heatmap(group_7_11_st_adata2, degs_7_11_tf, groupby='leiden', figsize=(5,7),
              #var_group_positions=[(0,10), (11, 12)], 
              vmin=-4, 
              vmax=4,
              cmap='vlag',
              var_group_rotation=0, dendrogram=False,
             #standard_scale ='var',
              swap_axes =True,
              show_gene_labels = False,
              save = 'heatmap_7_11_tf_scale.pdf'
             )

sc.tl.rank_genes_groups(st_adata2,groupby='leiden', groups = ['pMN_2'],reference ='pMN_4', method='t-test',n_genes =20000 )
x = sc.get.rank_genes_groups_df(st_adata2, group="pMN_2"); x.to_csv('pMN_4_pMN_2.txt',sep='\t')
degs_4_2 = get_DEG('pMN_4_pMN_2.txt')
cells = st_adata2.obs
group_4_2_cells = cells[(cells['leiden']=='pMN_4')|(cells['leiden']=='pMN_2')].index
group_4_2_st_adata2 = st_adata2[group_4_2_cells,:]
sc.pp.scale(group_4_2_st_adata2)
degs_4_2 = [g.encode('utf-8') for g in degs_4_2]
sc.pl.heatmap(group_4_2_st_adata2, degs_4_2, groupby='leiden', figsize=(5,7),
              #var_group_positions=[(0,10), (11, 12)], 
              vmin=-2,
              vmax=2,
              cmap='vlag',
              var_group_rotation=0, dendrogram=False,
             #standard_scale ='var',
              swap_axes =True,
              show_gene_labels = False,
              save = 'heatmap_4_2_scale.pdf'
             )
degs_4_2_tf = get_tf('pMN_4_pMN_2.txt',degs_4_2)
sc.pl.heatmap(group_4_2_st_adata2, degs_4_2_tf, groupby='leiden', figsize=(5,8),
              #var_group_positions=[(0,10), (11, 12)], 
              vmin=-3, 
              vmax=4,
              cmap='vlag',
              var_group_rotation=0, dendrogram=False,
             #standard_scale ='var',
              swap_axes =True,
              show_gene_labels = True,
              save = 'heatmap_4_2_tf_scale.pdf'
             )


sc.pl.stacked_violin(group_4_2_st_adata2, degs_4_2_tf, groupby='leiden',
                         var_group_positions=[(7, 8)], swap_axes =True, figsize=(3,40),
                     save='violin_4_2_tf.pdf'
                    )
sc.settings.set_figure_params(dpi=100,dpi_save=100,frameon=True)
sc.tl.rank_genes_groups(st_adata2,groupby='leiden', groups = ['pMN > Motor neuron_7'],reference ='pMN_2', method='t-test',n_genes =20000 )

x = sc.get.rank_genes_groups_df(st_adata2, group="pMN > Motor neuron_7"); x.to_csv('pMN2_pMNMotor7.txt',sep='\t')

degs_2_7 = get_DEG('pMN2_pMNMotor7.txt')



cells = st_adata2.obs
group_2_7_cells = cells[(cells['leiden']=='pMN > Motor neuron_7')|(cells['leiden']=='pMN_2')].index
group_2_7_st_adata2 = st_adata2[group_2_7_cells,:]



degs_2_7 = [g.encode('utf-8') for g in degs_2_7]



sc.pp.scale(group_2_7_st_adata2)



sc.pl.heatmap(group_2_7_st_adata2, degs_2_7, groupby='leiden', figsize=(5,7),
              #var_group_positions=[(0,10), (11, 12)], 
              vmin=-2, 
              vmax=2,
              cmap='vlag',
              var_group_rotation=0, dendrogram=False,
              #standard_scale ='var',
              swap_axes =True,
              show_gene_labels = False,
              save = 'heatmap_2_7_scale.pdf'
             )




degs_2_7_tf = get_tf('0721\pMN2_pMNMotor7.txt',degs_2_7)




sc.pl.heatmap(group_2_7_st_adata2,degs_2_7_tf, groupby='leiden', figsize=(5,8),
              #var_group_positions=[(0,10), (11, 12)], 
              vmin=-4, 
              vmax=4,
              cmap='vlag',
              var_group_rotation=0, dendrogram=False,
             #standard_scale ='var',
              swap_axes =True,
              show_gene_labels = False,
              save = 'heatmap_2_7_tf_scale.pdf'
             )
print(len(degs_2_7_tf))


###############################



sc.settings.set_figure_params(dpi=100,dpi_save=1000,frameon=True)
sc.tl.rank_genes_groups(st_adata2,groupby='leiden', groups = ['OPC_11'],reference ='pMN_2', method='t-test',n_genes =20000 )



x = sc.get.rank_genes_groups_df(st_adata2, group="OPC_11"); x.to_csv('pMN2_OPC11.txt',sep='\t')



degs_2_11 = get_DEG('pMN2_OPC11.txt')



cells = st_adata2.obs
group_2_11_cells = cells[(cells['leiden']=='OPC_11')|(cells['leiden']=='pMN_2')].index
group_2_11_st_adata2 = st_adata2[group_2_11_cells,:]



sc.pp.scale(group_2_11_st_adata2)



degs_2_11 = [g.encode('utf-8') for g in degs_2_11]



sc.pl.heatmap(group_2_11_st_adata2, degs_2_11, groupby='leiden', figsize=(5,7),
              #var_group_positions=[(0,10), (11, 12)], 
              vmin=-2, 
              vmax=2,
              cmap='vlag',
              var_group_rotation=0, dendrogram=False,
              #standard_scale ='var',
              swap_axes =True,
              show_gene_labels = False,
              save = 'heatmap_2_11_scale.pdf'
             )



degs_2_11_tf = get_tf('pMN2_OPC11.txt',degs_2_11)



degs_2_11_tf = [g.encode('utf-8') for g in degs_2_11_tf]



sc.pl.heatmap(group_2_11_st_adata2,degs_2_11_tf, groupby='leiden', figsize=(5,8),
              #var_group_positions=[(0,10), (11, 12)], 
              vmin=-2, 
              vmax=2,
              cmap='vlag',
              var_group_rotation=0, dendrogram=False,
             #standard_scale ='var',
              swap_axes =True,
              show_gene_labels = False,
              save = 'heatmap_2_11_tf_scale.pdf'
             )


############################# curve plot




st_adata2.obs['dpt_pseudotime'] = adata3.obs['dpt_pseudotime']




st_adata2[['AAACCCAAGAAATTGC-1'],:].obs['leiden']




cells = st_adata2.obs
core_cells = cells[(cells['leiden']=='pMN > Motor neuron_7')|
                         (cells['leiden']=='pMN_4')|(cells['leiden']=='pMN_2')|
                        (cells['leiden']=='Motor neuron_1')|(cells['leiden']=='Motor neuron_0')|(cells['leiden']=='Motor neuron_5')].index
core_st_adata2 = st_adata2[core_cells,:]




def get_sort_DEG(f):
    #f = 'pMNMotor7_OPC11.txt'
    df = pd.read_table(f,index_col=0);df.index = df['names']
    df = df.sort_values(by='pvals_adj')
    ups = df[(df['pvals_adj']<0.001)&(df['logfoldchanges']>1)].index.tolist()
    print(len(ups))
    downs = df[(df['pvals_adj']<0.001)&(df['logfoldchanges']<-1)].index.tolist()
    print(len(downs))
    np.random.shuffle(ups);np.random.shuffle(downs)
    #degs = ups + downs
    #fo = f.replace('.txt','_DEGs.txt')
    #with open(fo,'w') as o:
    #    for gene in degs:
    #        o.write(gene+'\n')
    return ups,downs




def get_plot_dict(df):
    #df = df.iloc[:,8:]
    df = df.sort_values(by='pseudotime')
    xy_dict = {}
    genes = df.columns[:-1]
    for gene in genes:
        ys = df[gene]
        ys = ys[ys>0] ###
        xs = df['pseudotime'][ys.index]
        #print(np.mean(ys))
        ysmoothed = gaussian_filter1d(ys.tolist(), sigma=80)
        #print(max(ysmoothed))
        xy_dict[gene] = [xs,ysmoothed]
    return xy_dict




branch_ups,branch_downs = get_sort_DEG('0721\pMN2_pMNMotor7.txt')


up_df = core_st_adata2[:,branch_ups[:10]].to_df()
up_df['pseudotime'] = core_st_adata2.obs['dpt_pseudotime']
xy_dict = get_plot_dict(up_df)




fig = plt.figure(figsize=(12,8))
colors=plt.get_cmap('tab10').colors[:10]
for gene,color in zip(xy_dict.keys(),colors):
    xs,ysmoothed = xy_dict[gene]
    #print(max(ysmoothed))
    plt.plot(xs, ysmoothed,label=gene,color=color,linewidth=3)

plt.grid(False)
sns.despine(right=True,top=True)
plt.xlabel('Pseudotime')
plt.ylabel('log2(TPM+1)')
plt.legend(loc='best')
plt.legend(frameon=False) 
plt.savefig('branch1_up_DEGs.pdf')




down_df = core_st_adata2[:,branch_downs[:10]].to_df()
down_df['pseudotime'] = core_st_adata2.obs['dpt_pseudotime']
xy_dict = get_plot_dict(down_df)




fig = plt.figure(figsize=(12,8))
colors=plt.get_cmap('tab10').colors[:10]
for gene,color in zip(xy_dict.keys(),colors):
    xs,ysmoothed = xy_dict[gene]
    #print(max(ysmoothed))
    plt.plot(xs, ysmoothed,label=gene,color=color,linewidth=3)

plt.grid(False)
sns.despine(right=True,top=True)
plt.xlabel('Pseudotime')
plt.ylabel('log2(TPM+1)')
plt.legend(loc='best')
plt.legend(frameon=False) 
plt.savefig('branch1_down_DEGs.pdf')




def get_lineplot_tf(degs):
    degs_tf = []
    for gene in degs:
        if gene in tfs:
            degs_tf.append(gene)
    return degs_tf



ups_tf = get_lineplot_tf(branch_ups)
up_df = core_st_adata2[:,ups_tf[:10]].to_df()
up_df['pseudotime'] = core_st_adata2.obs['dpt_pseudotime']
xy_dict = get_plot_dict(up_df)




fig = plt.figure(figsize=(12,8))
colors=plt.get_cmap('tab10').colors[:10]
for gene,color in zip(xy_dict.keys(),colors):
    xs,ysmoothed = xy_dict[gene]
    #print(max(ysmoothed))
    plt.plot(xs, ysmoothed,label=gene,color=color,linewidth=3)

plt.grid(False)
sns.despine(right=True,top=True)
plt.xlabel('Pseudotime')
plt.ylabel('log2(TPM+1)')
plt.legend(loc='best')
plt.legend(frameon=False) 
plt.savefig('branch1_up_TFs.pdf')




downs_tf = get_lineplot_tf(branch_downs)[:7]
downs_tf+=['tox3','myt1a','myt1b']
down_df = core_st_adata2[:,downs_tf[:10]].to_df()
down_df['pseudotime'] = core_st_adata2.obs['dpt_pseudotime']
xy_dict = get_plot_dict(down_df)




fig = plt.figure(figsize=(12,8))
colors=plt.get_cmap('tab10').colors[:10]
for gene,color in zip(xy_dict.keys(),colors):
    xs,ysmoothed = xy_dict[gene]
    #print(max(ysmoothed))
    plt.plot(xs, ysmoothed,label=gene,color=color,linewidth=3)

plt.grid(False)
sns.despine(right=True,top=True)
plt.xlabel('Pseudotime')
plt.ylabel('log2(TPM+1)')
plt.legend(loc='best')
plt.legend(frameon=False) 
plt.savefig('branch1_down_TFs_3.pdf')




sc.pl.stacked_violin(core_st_adata2, ['lbr','myt1a','myt1b','tox3'], groupby='leiden',
                         var_group_positions=[(7, 8)], swap_axes =True, figsize=(8,4),
                     save='branch1_violin.pdf'
                    )




cells = st_adata2.obs
core_cells = cells[(cells['leiden']=='pMN_4')|(cells['leiden']=='pMN_2')|(cells['leiden']=='OPC_11')].index
core_st_adata2 = st_adata2[core_cells,:]



def get_sort_DEG(f):
    #f = 'pMNMotor7_OPC11.txt'
    df = pd.read_table(f,index_col=0);df.index = df['names']
    df = df.sort_values(by='pvals_adj')
    ups = df[(df['pvals_adj']<0.001)&(df['logfoldchanges']>1)].index.tolist()
    print(len(ups))
    downs = df[(df['pvals_adj']<0.001)&(df['logfoldchanges']<-1)].index.tolist()
    print(len(downs))
    np.random.shuffle(ups);np.random.shuffle(downs)
    #degs = ups + downs
    #fo = f.replace('.txt','_DEGs.txt')
    #with open(fo,'w') as o:
    #    for gene in degs:
    #        o.write(gene+'\n')
    return ups,downs


def get_plot_dict(df):
    #df = df.iloc[:,8:]
    df = df.sort_values(by='pseudotime')
    xy_dict = {}
    genes = df.columns[:-1]
    for gene in genes:
        ys = df[gene]
        ys = ys[ys>0] ###
        xs = df['pseudotime'][ys.index]
        #print(np.mean(ys))
        ysmoothed = gaussian_filter1d(ys.tolist(), sigma=80)
        #print(max(ysmoothed))
        xy_dict[gene] = [xs,ysmoothed]
    return xy_dict


branch_ups,branch_downs = get_sort_DEG('0721\pMN2_OPC11.txt')




up_df = core_st_adata2[:,branch_ups[:10]].to_df()
up_df['pseudotime'] = core_st_adata2.obs['dpt_pseudotime']
xy_dict = get_plot_dict(up_df)




fig = plt.figure(figsize=(12,8))
colors=plt.get_cmap('tab10').colors[:10]
for gene,color in zip(xy_dict.keys(),colors):
    xs,ysmoothed = xy_dict[gene]
    #print(max(ysmoothed))
    plt.plot(xs, ysmoothed,label=gene,color=color,linewidth=3)

plt.grid(False)
sns.despine(right=True,top=True)
plt.xlabel('Pseudotime')
plt.ylabel('log2(TPM+1)')
plt.legend(loc='best')
plt.legend(frameon=False) 
plt.savefig('branch2_up_DEGs.pdf')




down_df = core_st_adata2[:,branch_downs[:10]].to_df()
down_df['pseudotime'] = core_st_adata2.obs['dpt_pseudotime']
xy_dict = get_plot_dict(down_df)



fig = plt.figure(figsize=(12,8))
colors=plt.get_cmap('tab10').colors[:10]
for gene,color in zip(xy_dict.keys(),colors):
    xs,ysmoothed = xy_dict[gene]
    #print(max(ysmoothed))
    plt.plot(xs, ysmoothed,label=gene,color=color,linewidth=3)

plt.grid(False)
sns.despine(right=True,top=True)
plt.xlabel('Pseudotime')
plt.ylabel('log2(TPM+1)')
plt.legend(loc='best')
plt.legend(frameon=False) 
plt.savefig('branch2_down_DEGs.pdf')



def get_lineplot_tf(degs):
    degs_tf = []
    for gene in degs:
        if gene in tfs:
            degs_tf.append(gene)
    return degs_tf


ups_tf = get_lineplot_tf(branch_ups)
up_df = core_st_adata2[:,ups_tf[:10]].to_df()
up_df['pseudotime'] = core_st_adata2.obs['dpt_pseudotime']
xy_dict = get_plot_dict(up_df)




fig = plt.figure(figsize=(12,8))
colors=plt.get_cmap('tab10').colors[:10]
for gene,color in zip(xy_dict.keys(),colors):
    xs,ysmoothed = xy_dict[gene]
    #print(max(ysmoothed))
    plt.plot(xs, ysmoothed,label=gene,color=color,linewidth=3)

plt.grid(False)
sns.despine(right=True,top=True)
plt.xlabel('Pseudotime')
plt.ylabel('log2(TPM+1)')
plt.legend(loc='best')
plt.legend(frameon=False) 
plt.savefig('branch2_up_TFs.pdf')




downs_tf = get_lineplot_tf(branch_downs)[:7]
downs_tf+=['tox3','myt1a','myt1b']
down_df = core_st_adata2[:,downs_tf[:10]].to_df()
down_df['pseudotime'] = core_st_adata2.obs['dpt_pseudotime']
xy_dict = get_plot_dict(down_df)



fig = plt.figure(figsize=(12,8))
colors=plt.get_cmap('tab10').colors[:10]
for gene,color in zip(xy_dict.keys(),colors):
    xs,ysmoothed = xy_dict[gene]
    #print(max(ysmoothed))
    plt.plot(xs, ysmoothed,label=gene,color=color,linewidth=3)

plt.grid(False)
sns.despine(right=True,top=True)
plt.xlabel('Pseudotime')
plt.ylabel('log2(TPM+1)')
plt.legend(loc='best')
plt.legend(frameon=False) 
plt.savefig('branch2_down_TFs_3.pdf')




sc.pl.stacked_violin(core_st_adata2, ['lbr','myt1a','myt1b','tox3'], groupby='leiden',
                         var_group_positions=[(7, 8)], swap_axes =True, figsize=(4,4),
                     save='branch2_violin.pdf'
                    )