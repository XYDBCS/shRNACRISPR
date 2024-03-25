import pandas as pd
import matplotlib.pyplot as plt
# import seaborn as sns
import numpy as np
import seaborn as sns
#############################################################
'''
Draw the high expression essentail genes and low expression essential gens
'''
##############################################################

'''
filename0 = './data/Gene Expression 530 Cell lines.xlsx'
Exp = pd.read_excel(filename0, index_col=0,engine='openpyxl')
Exp.sort_values('GeneID', inplace=True)
genes0 = Exp.index.tolist()   #genes with expression info

filename1 = "./data/Essentiality CRISPR 530 Cell lines.xlsx"

Crispr = pd.read_excel(filename1, index_col=0,engine='openpyxl')
genes1 = Crispr.index.tolist()  #genes with crispr info

# filename2 = "./data/Essentiality shRNA 530 Cell lines.xlsx"
# shRNA = pd.read_excel(filename2, index_col=0,engine='openpyxl')
#
# n = shRNA.shape[1]
# for i in range(n):
#     tep_row = shRNA.iloc[:,i].values
#     if np.count_nonzero(tep_row>0)+np.count_nonzero(tep_row<0) <3000:
#         shRNA.drop(1,i)
# print("shRNA", shRNA.shape[0])

filename2 = "./data/Essentiality shRNA 530 Cell lines.xlsx"

shRNA = pd.read_excel(filename2, index_col=0,engine='openpyxl')
shRNA = shRNA.fillna(0)
# print("shRNA",shRNA)

n0_shRNA = (shRNA != 0).astype(int).sum(axis=0)
shRNA3000 = n0_shRNA[n0_shRNA.astype(float) > 3000]
# print(shRNA3000)

cells = shRNA3000.index.tolist()
print("cells",cells)
shRNA = shRNA.reindex(cells, axis='columns')
# print(shRNA)


genes2 = shRNA.index.tolist()   #genes with shRNA info
genes = list(set(genes1).intersection(set(genes2)).intersection(set(genes0)))


Exp_re = Exp.loc[genes].drop(columns=['GeneID'])
# Exp_re['exp_sum'] = Exp_re.apply(lambda x: x.sum(), axis=1) #按行求和，添加为新列
# Exp_re.sort_values("exp_sum", inplace=True)
genes = Exp_re.index.tolist()  #update order of genes accroding to expression


# cells = Crispr.columns.tolist()
# del(cells[0])
#print(cells)

Crispr_re = Crispr.loc[genes].drop(columns=['GeneID'])
shRNA_re = shRNA.loc[genes].drop(columns=['GeneID']).fillna(0)

#uniform header(cell line order)of 3 forms
Exp_re = Exp_re.reindex(cells, axis = 'columns')
shRNA_re = shRNA_re.reindex(cells, axis = 'columns').fillna(0)
Crispr_re = Crispr_re.reindex(cells, axis = 'columns').fillna(0)
print(Exp_re, Crispr_re, shRNA_re)
'''
###############################calculate the correlation values###########################
###task 3#################################################################################

filename1 = "./data/Essentiality CRISPR 254 cell lines.xlsx"
Crispr = pd.read_excel(filename1, index_col=0)
cell_lines = Crispr.columns.tolist()
Crispr = Crispr.reindex(cell_lines, axis='columns')
gene_crispr = Crispr.index.tolist()

filename2 = "./data/Essentiality shRNA 254 cell lines.xlsx"
shRNA = pd.read_excel(filename2, index_col=0)
shRNA = shRNA.fillna(0)
shRNA = shRNA.reindex(cell_lines, axis='columns')

filename0 = './data/Gene Expression 530 Cell lines.xlsx'
Exp = pd.read_excel(filename0, index_col=0)
genes0 = Exp.index.tolist()   #genes with expression info
Exp = Exp.reindex(cell_lines, axis = 'columns')
gene_all = Exp.index.tolist()
gene_shRNA = shRNA.index.tolist()

genes = list(set(gene_shRNA).intersection(set(gene_crispr)).intersection(set(gene_all)))
print(len(genes))
shRNA = shRNA.loc[genes].drop(columns=['GeneID'])
Crispr = Crispr.loc[genes].drop(columns=['GeneID'])
Exp = Exp.loc[genes].drop(columns=['GeneID'])

shRNA_v = shRNA.values
Crispr_v = Crispr.values
Exp_v = Exp.values


pearson_corr = []
pearson_kendall = []
pearson_spearman = []
cell_names = []
print("cell_n", )
for i in range(254):
    shRNA_df = pd.Series(shRNA_v[:,i])
    Crispr_df = pd.Series(Crispr_v[:,i])
    Exp_df = pd.Series(Exp_v[:,i])
    pearson_corr.append(shRNA_df.corr(Exp_df))
    pearson_kendall.append(Crispr_df.corr(Exp_df))
    pearson_spearman.append(shRNA_df.corr(Crispr_df))

#     l1,l2,l3 = (list(t)for t in zip(*sorted(zip(pearson_corr,pearson_kendall,pearson_spearman),reverse = True)))
#
#     # pearson_corr.append(np.corrcoef(shRNA_v[:,i], Exp_v[:,i])[0,1])
#     # pearson_kendall.append(np.corrcoef(Crispr_v[:,i],Exp_v[:,i])[0,1])
#     # pearson_spearman.append(np.corrcoef(shRNA_v[:,i], Crispr_v[:,i])[0,1])
#     # pearson_kendall.append(shRNA.iloc[:, [i]].corr(Crispr.iloc[:, [i]], method='kendall'))
#     # pearson_spearman.append(shRNA.iloc[:, [i]].corr(Crispr.iloc[:, [i]], method='spearman'))
#
# y = range(254)
# plt.plot(y,l1,  lw=1,color = 'darkviolet', label='shRNA_expression')
# plt.plot(y,l2,  lw=1, color = 'teal',  label='CRISPR_expression')
# plt.plot( y,l3, lw=1,color = 'red', label='shRNA_CRISPR')
#
# print("shRNA", np.average(np.array(pearson_corr)), "CRISPR", np.average(np.array(pearson_kendall)), "shRNA_CRISPR", np.average(np.array(pearson_spearman)))
#
#
#
# plt.title('the correlation between shRNA,CRISPR, and gene expression')
# plt.xlabel('cell lines')
# plt.ylabel('correlation value')
# plt.legend()
# plt.show()

#########################################################################################################
###task 3################################################################################################
###############################calculate the correlation values##########################################@


#######################################task 1: draw the distributation of genes and task 2: distributation of gene expression task3: correlation###################################
####################################################################################################
top_rate = 0.3#top expression

top_rate_pre = 0.3#depend on the prediction essential genes

top_n = int(top_rate_pre*15029)
top_exp_rate = top_rate
filename1 = "./data/Essentiality CRISPR 254 cell lines.xlsx"
# filename11 = "./Unsupervised_prediction score_RQ_0.1_0.0001_256_6000.txt"
# filename11 ="./unsupervised_prediction score_0.15_0.0001_norm_6000_256.txt"
filename11 ="./final/Plato_local/Unsupervised_1000_CRISPR.txt"
# Crispr = pd.read_excel(filename1, index_col=0,engine='openpyxl')
Crispr = pd.read_excel(filename1, index_col=0)
cell_lines = Crispr.columns.tolist()
pred_Crispr = np.loadtxt(filename11)
Crispr = Crispr.reindex(cell_lines, axis='columns')
Crispr['score_ess'] = pred_Crispr[:,0]
Crispr['score_non_ess'] = pred_Crispr[:,2]
gene_crispr = Crispr.index.tolist()

filename2 = "./data/Essentiality shRNA 254 cell lines.xlsx"
# filename22 = "./shRNA_unsupervised_Del_prediction score_0.1_0.0001_norm_256_6000.txt"
filename22 ="./final/Plato_local/shRNA_Unsupervised_1000.txt"

# shRNA = pd.read_excel(filename2, index_col=0,engine='openpyxl')
shRNA = pd.read_excel(filename2, index_col=0)
shRNA = shRNA.fillna(0)
pred_shRNA = np.loadtxt(filename22)
shRNA = shRNA.reindex(cell_lines, axis='columns')
shRNA['score_ess'] = pred_shRNA[:,0]
shRNA['score_non_ess'] = pred_shRNA[:,2]



filename0 = './data/Gene Expression 530 Cell lines.xlsx'
# Exp = pd.read_excel(filename0, index_col=0,engine='openpyxl')
Exp = pd.read_excel(filename0, index_col=0)
genes0 = Exp.index.tolist()   #genes with expression info
Exp = Exp.reindex(cell_lines, axis = 'columns')
gene_all = Exp.index.tolist()
gene_shRNA = shRNA.index.tolist()

genes = list(set(gene_shRNA).intersection(set(gene_crispr)).intersection(set(gene_all)))
print(len(genes))
shRNA = shRNA.loc[genes].drop(columns=['GeneID'])
Crispr = Crispr.loc[genes].drop(columns=['GeneID'])
Exp = Exp.loc[genes].drop(columns=['GeneID'])

#########################################
Crispr.sort_values("score_ess",  ascending=False, inplace=True)
gene_C_n = Crispr.shape[0]
print("gene_#", gene_C_n, Crispr.shape)
label = np.zeros([gene_C_n])
label[0:top_n] = 1
Crispr['label_ess'] = label
Crispr.sort_values("score_non_ess",  ascending=False, inplace=True)
Crispr['label_ness'] = label
print("sum,",sum(Crispr['label_ess'][0:2000]))
print(Crispr)
Crispr = Crispr.loc[genes]


shRNA.sort_values('score_ess',  ascending=False, inplace = True)
gene_s_n = shRNA.shape[0]
print("gene#_shRNA", gene_s_n)
label_s = np.zeros([gene_s_n])
label_s[0:top_n] = 1
shRNA['label_ess'] = label_s
shRNA.sort_values('score_non_ess', ascending=False, inplace = True)
shRNA['label_ness'] = label_s
shRNA = shRNA.loc[genes]
print(shRNA)
#print("sum,",sum(shRNA['label_ess'][0:2000]))



exp_num = Exp.values
Crispr_num = Crispr.values
print("Crispr_num", Crispr_num)
shRNA_num = shRNA.values
m,n = exp_num.shape

CRISP_e = np.zeros([n])
CRISP_ne = np.zeros([n])
CRISP_e_S = np.zeros([n])
CRISP_ne_S = np.zeros([n])

shRNA_e = np.zeros([n])
shRNA_ne = np.zeros([n])
shRNA_e_S = np.zeros([n])
shRNA_ne_S = np.zeros([n])


#a is the numpy array, when k==0, return all the >0 numbre and index, when k==1, return all the <0 numbers and index
# def index_larger_small (a,k):
#     d = np.arange(len(a))
#     b_ne = d[a>0]
#     c_e = d[a<0]
#     if k ==0:#return essential
#         return a[c_e], c_e
#     else: #return non_essential
#         return a[b_ne], b_ne
#
# print("m,n", m, n, Crispr_num.shape)
# # print("Crispr_num_n+2", np.sum(Crispr_num[:, n+3]),  np.sum(Crispr_num[:, n+2]))
corr_non_shRNA = []
corr_shRNA = []
corre_Crispr = []
corrn_Crispr = []

Crispr_essS_exp = []
shRNA_essS_exp = []
for i in range(n):
    sort_exp = np.argsort(exp_num[:,i])
    exp_top = int((top_exp_rate) * m)
    gene_exp_25 = sort_exp[(m-exp_top):]
    gene_exp_25 = sort_exp[:]
    CRISP_e_S[i] = np.sum(Crispr_num[:,i]<0)
    CRISP_ne_S[i] = np.sum(Crispr_num[:,i]>0)
    CRISP_e[i] = np.sum(Crispr_num[gene_exp_25,i]<0)/CRISP_e_S[i]
    CRISP_ne[i] = np.sum(Crispr_num[gene_exp_25, i] > 0)/CRISP_ne_S[i]

# #  ##########################Calculte the expression of all essential genes in CRISPR and shRNA ######################
#
#     correlation_Crisp_essS_exp = []
#     correlation_Crisp_nssS_exp = []
#     # CRISP_i_sort = np.argsort(Crispr_num[:, i])
#     CRISP_i_sort = np.argsort(-Crispr_num[:, i])
#     th_e = Crispr_num[CRISP_i_sort[int(top_rate * m)], i]
#     for j in range(m):
#         if  Crispr_num[j, i] < th_e and Crispr_num[j, n+2]==1:
#         # if Crispr_num[j, i] > th_e and Crispr_num[j, n + 3] == 1:
#             Crispr_essS_exp.append([Crispr_num[j, i], exp_num[j, i]])
#             correlation_Crisp_essS_exp.append([Crispr_num[j, i], exp_num[j, i]])
#         elif Crispr_num[j, i] > th_e and Crispr_num[j, n + 3] == 1:
#             correlation_Crisp_nssS_exp.append([Crispr_num[j, i], exp_num[j, i]])
#
#
#
#     correlation_Crisp_essS_exp = np.array(correlation_Crisp_essS_exp)
#     correlation_Crisp_nssS_exp = np.array(correlation_Crisp_nssS_exp)
#     corrn_Crispr.append(np.corrcoef(correlation_Crisp_nssS_exp[:,0], correlation_Crisp_nssS_exp[:,1])[0,1])
#     corre_Crispr.append(np.corrcoef(correlation_Crisp_essS_exp[:,0], correlation_Crisp_essS_exp[:,1])[0,1])
#
#
#
#
#     correlation_shRNA_essS_exp = []
#     correlation_shRNA_nssS_exp = []
#
#     # shRNA_i_sort = np.argsort(shRNA_num[:, i])
#     shRNA_i_sort = np.argsort(-shRNA_num[:, i])
#     th_S_e = shRNA_num[shRNA_i_sort[int(top_rate * m)], i]
#     for j in range(m):
#          if shRNA_num[j,i] < th_S_e and shRNA_num[j, n+2] ==1:
#          # if shRNA_num[j, i] >th_S_e and shRNA_num[j, n + 3] == 1:
#             shRNA_essS_exp.append([shRNA_num[j,i], exp_num[j,i]])
#             correlation_shRNA_essS_exp.append([shRNA_num[j,i], exp_num[j,i]])
#          elif shRNA_num[j, i] >th_S_e and shRNA_num[j, n + 3] == 1:
#              correlation_shRNA_nssS_exp.append([shRNA_num[j,i], exp_num[j,i]])
#
#     correlation_shRNA_nssS_exp = np.array(correlation_shRNA_nssS_exp)
#     correlation_shRNA_essS_exp = np.array(correlation_shRNA_essS_exp)
#     corr_shRNA.append(np.corrcoef(correlation_shRNA_essS_exp[:,0], correlation_shRNA_essS_exp[:,1])[0,1])
#     corr_non_shRNA.append(np.corrcoef(correlation_shRNA_nssS_exp[:,0], correlation_shRNA_nssS_exp[:,1])[0,1])
# print("avg_shRNA_ess", np.average(np.array(corr_shRNA)), np.average(np.array(corr_non_shRNA)))
# print("avg_CRISPR_ess", np.average(np.array(corre_Crispr)), np.average(np.array(corrn_Crispr)))
#
# l1,l2,l3,l4 = (list(t) for t in zip(*sorted(zip(corr_shRNA, corr_non_shRNA, corre_Crispr, corrn_Crispr), reverse=True)))
# plt.plot(range(n),   l1,  lw=1,color = 'darkviolet', label='essential: shRNA_expression')
# plt.plot(range(n), l2,  lw=1,linestyle='--', color = 'darkviolet',  label=' nonessential: shRNA_expression')#teal
# plt.plot(range(n), l3,  lw=1, color = 'teal',  label='essential: CRISPR_expression')#darkviolet
# plt.plot(range(n), l4,  lw=1,linestyle='--', color = 'teal',  label='nonessential: CRISPR_expression')
#
# plt.title('the correlation between essential (nonessential) scores and gene expression')
# plt.xlabel('cell lines')
# plt.ylabel('correlation value')
# plt.legend()
# plt.show()


##########################################correlation##########################################################

#
# Crispr_essS_exp = np.array(Crispr_essS_exp)
# shRNA_essS_exp = np.array(shRNA_essS_exp)
# print("Crispr_list", Crispr_essS_exp.shape, Crispr_essS_exp)
# print("shRNA_list", shRNA_essS_exp.shape, shRNA_essS_exp)
#
#
# df_shRNA = pd.DataFrame(shRNA_essS_exp,index=None,columns = ['score','expression'])
# df_Crispr = pd.DataFrame(Crispr_essS_exp, index=None, columns=['score', 'expression'])
# df_shRNA['label'] = 'shRNA'
# df_Crispr['label'] = 'Crispr'
# # sns.jointplot(x='score',y='expression',data=df_Crispr,hue='label')
# # sns.jointplot(x='score',y='expression',data=df_shRNA,hue='label')
#
# data = pd.concat([df_shRNA, df_Crispr],sort=False,axis=0,ignore_index=True)
# data = data.loc[~(data==0).all(axis=1)]
# sns.jointplot(x='score',y='expression',data=data,hue='label',  s = 35, palette={
#     # 'shRNA': 'blueviolet',
#     # 'Crispr': 'firebrick'
#     'shRNA': 'blueviolet',
#     'Crispr': 'crimson'
# })
#
#
#
#
#
# # plt.scatter(Crispr_essS_exp[:,0], Crispr_essS_exp[:,1], s=6, label="Cripsr", alpha=0.6)
# # plt.scatter(shRNA_essS_exp[:,0], shRNA_essS_exp[:,1], s=6, label="shRNA", alpha=0.6)
# #
# # plt.xlabel("essential score")
# # plt.ylabel("expression")
# # plt.legend()
# # plt.title("average_essential genes")
# plt.show()

 ##########################Calculte the expression of all essential genes in CRISPR and shRNA ######################

#Draw the high expression and low expression figures    
###########CRISPR#####################################
    CRISP_i_sort = np.argsort(Crispr_num[:,i])
    th_e = Crispr_num[CRISP_i_sort[int(top_rate*m)],i]
    # CRISP_e[i] = np.sum(Crispr_num[gene_exp_25, i] < th_e)
    for j in range(exp_top):
        if Crispr_num[gene_exp_25[j], i] < th_e and Crispr_num[gene_exp_25[j], n+2]==1:
            CRISP_e[i] = CRISP_e[i] + 1
    th_ne = Crispr_num[CRISP_i_sort[int((1-top_rate)*m)],i]
    # CRISP_ne[i] = np.sum(Crispr_num[gene_exp_25, i] > th_ne)
    for j in range(exp_top):
        if Crispr_num[gene_exp_25[j],i]>th_ne and Crispr_num[gene_exp_25[j], n+3]==1:
            CRISP_ne[i] = CRISP_ne[i]+1
#####shRNA###########################################
    shRNA_i_sort = np.argsort(shRNA_num[:, i])
    th_S_e = shRNA_num[shRNA_i_sort[int(top_rate * m)], i]
    th_S_ne = shRNA_num[shRNA_i_sort[int((1 - top_rate) * m)], i]

    # shRNA_e[i] = np.sum(shRNA_num[gene_exp_25, i] < th_S_e)
    
    for j in range(exp_top):
        if shRNA_num[gene_exp_25[j], i] < th_S_e and shRNA_num[gene_exp_25[j], n + 2] == 1:
            shRNA_e[i] = shRNA_e[i] + 1
           
    for j in range(exp_top):
        if shRNA_num[gene_exp_25[j], i] > th_S_ne and shRNA_num[gene_exp_25[j], n + 3] == 1:
            shRNA_ne[i] = shRNA_ne[i] + 1
def plot_line(y1, y3, y11, y33):
    print("y1.shape", y1.shape[0])
    x = range(y1.shape[0])
    plt.plot(x, y1, lw=1,color = 'teal', label='CRISPR_essential')
    plt.plot(x, y11, lw=1, color = 'teal', linestyle='dashed', label='CRISPR_non_essential')
    print("y11", y11)
    plt.plot(x, y3, lw=1,color = 'indigo', label='shRNA_essential')
    plt.plot(x, y33, lw=1,  color = 'darkviolet', linestyle='dashed', label='shRNA_non_essential')
    plt.title('top 10% high expression genes')
    # plt.title(' The number of essential genes and non essentical genes in different cell_lines')
    plt.xlabel('cell lines')
    plt.ylabel('essential #')
    #plt.xticks(x, genes, rotation=45)
    plt.legend()
    plt.show()
y1 = CRISP_e
y11 = CRISP_ne
y3 = shRNA_e
y33 = shRNA_ne

y3_ind = np.argsort(y3)
y1 = y1[y3_ind]
y3 = y3[y3_ind]
y11 = y11[y3_ind]
y33 = y33[y3_ind]

plot_line( y1, y3, y11, y33)

    # ess_CRISP_inter = np.intersect1d(Crispr_sort_cell[0:int(top_rate*m)], gene_exp_25)
    # print(ess_CRISP_inter)
    # CRISP_e[i] = np.sum(Crispr_num[ess_CRISP_inter, n+2] ==1)

    # CRISP_ne[i] = np.sum(Crispr_num[noness_CRISPR_inter, n+3] ==1)

    # ess_value, ess_index = index_larger_small(Crispr_num[gene_exp_25, i],0)
    # non_ess_value, non_index = index_larger_small(Crispr_num[gene_exp_25, i],1)
    # CRISP_e[i] = np.corrcoef(ess_value, exp_num[gene_exp_25[ess_index],i])[0,1]
    # CRISP_ne[i] = np.corrcoef(non_ess_value, exp_num[gene_exp_25[non_index],i]) [0,1]

    # shRNA_e_S[i] = np.sum(shRNA_num[:,i]<0)
    # shRNA_ne_S[i] = np.sum(shRNA_num[:,i]>0)

    # sh_ess_value, sh_ess_index = index_larger_small(shRNA_num[gene_exp_25,i], 0)
    # sh_non_ess_value, sh_non_index = index_larger_small(shRNA_num[gene_exp_25, i], 1)
    # shRNA_e[i] = np.sum(shRNA_num[gene_exp_25, i]<0)/shRNA_e_S[i]
    # shRNA_ne[i]=np.sum(shRNA_num[gene_exp_25, i]>0)/shRNA_ne_S[i]

    # shRNA_ne[i] = np.sum(shRNA_num[gene_exp_25, i] > th_S_ne)

    # shRNA_e[i] = np.sum(shRNA_num[ess_shRNA_inter, n+2] ==1)
    # shRNA_ne[i] = np.sum(shRNA_num[noness_shRNA_inter, n+3] ==1)


    # shRNA_e[i] = np.corrcoef(sh_ess_value, exp_num[gene_exp_25[sh_ess_index],i])[0,1]
    # shRNA_ne[i]= np.corrcoef(sh_non_ess_value, exp_num[gene_exp_25[sh_non_index],i])[0,1]


# Crispr_neg = (Crispr_re < 0).sum(axis=1)
# print('essential',Crispr_neg)
# Crispr_pos = (Crispr_re > 0).sum(axis=1)
# shRNA_neg = (shRNA_re < 0).sum(axis=1)
# shRNA_pos = (shRNA_re > 0).sum(axis=1)
#
#
# result = pd.concat([Crispr_neg, Crispr_pos, shRNA_neg, shRNA_pos], axis=1)
# result.columns = ['Crispr_essential','Crispr_non-essential', 'shRNA_essential','shRNA_non-essential']
#
# quarter = int(result.shape[0] * 0.95)
# top5 = result[quarter:]
# top5.sort_values('shRNA_essential', ascending=True, inplace=True)
#
# compare = top5[top5["shRNA_essential"] > top5["Crispr_essential"]]
# print('compare', compare)

# def plot_line(genes, y1, y3):
    # x = range(len(genes))
# def plot_line(y1, y3):
#     print("y1.shape", y1.shape[0])
#     x = range(y1.shape[0])
#     plt.plot(x, y1, lw=1, label='Crispr_essential')
#     plt.plot(x, y3, lw=1, label='shRNA_essential')
#     plt.title(' 25% low expression gene')
#     plt.xlabel('cell lines')
#     plt.ylabel('essential #')
#     #plt.xticks(x, genes, rotation=45)
#     plt.legend()
#     plt.show()

# names = top5.index.tolist()
# y1 = top5['Crispr_essential'].values.tolist()
# y3 = top5['shRNA_essential'].values.tolist()
# plot_line(names, y1, y3)

# names = top5.index.tolist()


# shRNA = np.argsort(shRNA_e_S)
# shRNA_e_S = shRNA_e_S[shRNA]
# shRNA_ne_S = shRNA_ne_S[shRNA]
# CRISP_e_S = CRISP_e_S[shRNA]
# CRISP_ne_S = CRISP_ne_S[shRNA]
# plot_line( CRISP_e_S, shRNA_e_S,CRISP_ne_S,shRNA_ne_S)