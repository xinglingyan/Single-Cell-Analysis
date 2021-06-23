# import packages

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

sc.settings.verbosity = 3         
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80)

# read file
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
#fig 1C
sc.pl.umap(adata3,color=['leiden','olig2','isl1'])
sc.pl.umap(adata3,color=['olig1',"olig2",'isl1','isl2a','lhx1a','sox10','nkx6.1', 'irx3a', 'LHX3','mnx2a','mnx2b','mnx1', 'gfap','fabp7a','vsx2','vsx1','gata3','tal1','gata2a','tal2','en1b'],save='markers.pdf')
sc.tl.rank_genes_groups(adata3, 'leiden', method='t-test',use_raw =False,n_genes =100)

marker_genes = {
    'pMN':{'olig1', 'olig2'},'Motor neuron':{'isl1', 'lhx1'},
    'OPC':{'olig1','olig2','sox10'}
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
#fig 1C
sc.pl.umap(adata3,color=['leiden'],save='Umap.pdf')

# trajectory
sc.tl.paga(adata3, groups='leiden',model='v1.2')
sc.settings.set_figure_params(dpi=100,dpi_save=1000,frameon=True)
#fig S1C
sc.pl.paga(adata3,color=['leiden'],save='abstract_paga.pdf')
sc.tl.draw_graph(adata3,init_pos='paga' ,layout ='fa')
sc.settings.set_figure_params(dpi=300,dpi_save=1000,frameon=True)
#fig S1A
sc.pl.draw_graph(adata3, color=['leiden'], legend_loc='on data', title='', legend_fontsize=6,size=4,frameon =True,save='.pdf')
#pseudotime
x = adata3.copy()
sc.tl.diffmap(x)
sc.pp.neighbors(x, n_neighbors=10)
sc.tl.draw_graph(x)
sc.settings.set_figure_params(dpi=300,dpi_save=1000,frameon=True)
sc.pl.draw_graph(x, color=['leiden','sox10'], legend_loc='on data', title='', legend_fontsize=6,size=4,frameon =True)
adata3.uns['iroot'] = np.flatnonzero(adata3.obs['leiden']  == 'pMN_4')[0]
sc.tl.dpt(adata3)#, n_dcs=5,min_group_size =0.05,n_branchings =1)
#fig S1B
sc.pl.draw_graph(adata3, color=['leiden','dpt_pseudotime'], legend_loc='on data', title='', 
    legend_fontsize=6,size=4,vmax = 0.5,color_map ='Blues',save='.pdf')
#save data
adata3.write('adata3.h5ad')
##########################################################################################
#fig S1E
cells = adata3.obs
G7_cells = cells[(cells['leiden']=='pMN_2')|(cells['leiden']=='pMN > Motor neuron_7')
|(cells['leiden']=='Motor neuron_1')|(cells['leiden']=='Motor neuron_0')|(cells['leiden']=='Motor neuron_5')].index
G7_all = adata3[G7_cells,:]
print(G7_all)
G7_cell_type = ['pMN_2','pMN > Motor neuron_7','Motor neuron_1','Motor neuron_0','Motor neuron_5']
sc.tl.rank_genes_groups(G7_all,groupby='leiden', groups = G7_cell_type ,reference ='rest', method='t-test',n_genes =20000 )
G7_genes = []
for cell_type in G7_cell_type:
    fo = 'G7_{}.txt'.format(cell_type)
    fo = fo.replace(' ',''); fo = fo.replace('>','')
    x = sc.get.rank_genes_groups_df(G7_all, group=cell_type);x.to_csv(fo,sep='\t')
    tmp_genes = x['names'][:100].tolist()
    for g in tmp_genes:
        if g in G7_genes:
            continue
        G7_genes.append(g)

sc.pp.scale(G7_all)
sc.pl.heatmap(G7_all, G7_genes, groupby='leiden', figsize=(8,12),
              #var_group_positions=[(0,10), (11, 12)], 
              vmin=-4, 
              vmax=4,
              cmap='vlag',
              var_group_rotation=0, dendrogram=False,
             #standard_scale ='var',
              swap_axes =True,
              show_gene_labels = False,
              save = 'heatmap.pdf'
             )
##########################################################################################
#fig 4
group_3_adata3 = adata3[group_3_cells,:]
sc.settings.set_figure_params(dpi=100,dpi_save=1000,frameon=True)
sc.tl.umap(group_3_adata3,min_dist=0.1,)
sc.pl.umap(group_3_adata3)
group_3_adata3.uns['leiden_colors']
sc.tl.leiden(group_3_adata3,resolution=0.24,directed=False)
sc.settings.set_figure_params(dpi=100,dpi_save=1000,frameon=True)
sc.pl.umap(group_3_adata3,color=['leiden','olig2','isl1','sox10'],save='.pdf')


##########################################################################################
#fig S2
sc.tl.tsne(adata3,)
sc.pl.tsne(adata3, color=['leiden'], legend_loc='on data', title='', legend_fontsize=6,size=4,frameon =True,save = '.pdf')
##########################################################################################
#fig S3
G5_cells = cells[(cells['leiden']=='pMN_2')|(cells['leiden']=='pMN_4')].index
G5_all = adata3[G5_cells,:]
G5_all.write('G5.h5ad')
G5_all.obs.to_csv('G5.index.txt',sep='\t')
G5_cell_type = ['pMN_2','pMN_4']
sc.tl.rank_genes_groups(G5_all,groupby='leiden', groups = G5_cell_type ,reference ='rest', method='t-test',n_genes =20000 )
G5_genes = []
for cell_type in G5_cell_type:
    fo = 'subgroupheatmap/G5_{}.txt'.format(cell_type)
    fo = fo.replace(' ',''); fo = fo.replace('>','')
    x = sc.get.rank_genes_groups_df(G5_all, group=cell_type);x.to_csv(fo,sep='\t')
    tmp_genes = x['names'][:100].tolist()
    for g in tmp_genes:
        if g in G5_genes:
            continue
        G5_genes.append(g)

sc.pp.scale(G5_all)
sc.pl.heatmap(G5_all, G5_genes, groupby='leiden', figsize=(8,12),
              #var_group_positions=[(0,10), (11, 12)], 
              vmin=-2, 
              vmax=2,
              #center=-0.5,
              cmap='vlag',
              var_group_rotation=0, dendrogram=False,
             #standard_scale ='var',
              swap_axes =True,
              #show_gene_labels = True,
              save = 'G5_heatmap.pdf'
             )
#
sc.pl.stacked_violin(G5_all, ['nhlh2','hmgb3b','hmga1a','ebf2','neurod1'], groupby='leiden',
                         var_group_positions=[(7, 8)], swap_axes =True, figsize=(3,2),save='myt1a_violin.pdf')

##########################################################################################

#fig S4
tmp_cells = cells[(cells['leiden']=='pMN > Motor neuron_7')|(cells['leiden']=='OPC_11')|
                 (cells['leiden']=='pMN_4')|(cells['leiden']=='pMN_2')].index
sub_adata_hv = sub_adata_hv[tmp_cells,:]
hv = sc.pp.highly_variable_genes(sub_adata,n_top_genes=4000,inplace=False)
sub_adata_hv = sub_adata[:, hv.highly_variable]
sub_adata_hv.uns['leiden_colors'] = ['#279e68','#aa40fc','#b5bd61','#98df8a']
sc.pl.paga(sub_adata_hv,color=['leiden'],threshold=0.1,
           save='subsample'
          )
sc.pl.draw_graph(sub_adata_hv,color=['leiden','dpt_pseudotime'], legend_loc='on data', title='', 
                 legend_fontsize=6,size=4,vmax=1,arrows=False,
                 save='susample'
                )

##########################################################################################
#fig S6
group_7_11_cells = cells[(cells['leiden']=='pMN > Motor neuron_7')|(cells['leiden']=='OPC_11')].index
group_7_11_st_adata2 = st_adata2[group_7_11_cells,:]
sc.pp.scale(group_7_11_st_adata2)
sc.pl.heatmap(group_7_11_st_adata2, degs_7_11, groupby='leiden', figsize=(5,7),
              vmin=-4, 
              vmax=4,
              cmap='vlag',
              var_group_rotation=0, dendrogram=False,
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
              vmin=-4, 
              vmax=4,
              cmap='vlag',
              var_group_rotation=0, dendrogram=False,
              swap_axes =True,
              show_gene_labels = False,
              save = 'heatmap_7_11_tf_scale.pdf'
             )









