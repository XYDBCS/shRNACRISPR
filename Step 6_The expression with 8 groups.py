import pandas as pd
import matplotlib.pyplot as plt
# import seaborn as sns
import numpy as np
import seaborn as sns

#######################################task 1: draw the distributation of genes and task 2: distributation of gene expression ###################################
####################################################################################################
top_rate = 0.1
top_rate_pre = 0.1
top_n = int(top_rate_pre * 15029)
top_exp_rate = top_rate
filename1 = "./data/Essentiality CRISPR 254 cell lines.xlsx"
# filename11 = "./Unsupervised_prediction score_RQ_0.1_0.0001_256_6000.txt"
filename11 = "./unsupervised_prediction score_0.15_0.0001_norm_6000_256.txt"
# filename11 ="./final/Plato_local/Unsupervised_1000_CRISPR.txt"
# Crispr = pd.read_excel(filename1, index_col=0,engine='openpyxl')
Crispr = pd.read_excel(filename1, index_col=0)
cell_lines = Crispr.columns.tolist()  # cell line name from CRISPR
# print("cell_lines", cell_lines)

pred_Crispr = np.loadtxt(filename11)
Crispr = Crispr.reindex(cell_lines, axis='columns')
Crispr['score_ess'] = pred_Crispr[:, 0]
Crispr['score_non_ess'] = pred_Crispr[:, 2]
gene_crispr = Crispr.index.tolist()  # gene name from CRISPR
# print("gene_name", gene_crispr)

filename2 = "./data/Essentiality shRNA 254 cell lines.xlsx"
filename22 = "./shRNA_unsupervised_Del_prediction score_0.1_0.0001_norm_256_6000.txt"
# filename22 ="./final/Plato_local/shRNA_Unsupervised_1000.txt"

# shRNA = pd.read_excel(filename2, index_col=0,engine='openpyxl')
shRNA = pd.read_excel(filename2, index_col=0)
shRNA = shRNA.fillna(0)
pred_shRNA = np.loadtxt(filename22)
shRNA = shRNA.reindex(cell_lines, axis='columns')
shRNA['score_ess'] = pred_shRNA[:, 0]
shRNA['score_non_ess'] = pred_shRNA[:, 2]

filename0 = './data/Gene Expression 530 Cell lines.xlsx'
# Exp = pd.read_excel(filename0, index_col=0,engine='openpyxl')
Exp = pd.read_excel(filename0, index_col=0)
genes0 = Exp.index.tolist()  # genes with expression info
Exp = Exp.reindex(cell_lines, axis='columns')
gene_all = Exp.index.tolist()  # gene from expression
gene_shRNA = shRNA.index.tolist()  # gene from shRNA

f = open("cell_line_name.txt",'w')
for i in cell_lines:
    f.write(str(i) + '\n')

genes = list(set(gene_shRNA).intersection(set(gene_crispr)).intersection(set(gene_all)))
print(len(genes))
shRNA = shRNA.loc[genes].drop(columns=['GeneID'])
Crispr = Crispr.loc[genes].drop(columns=['GeneID'])
Exp = Exp.loc[genes].drop(columns=['GeneID'])

#########################################
Crispr.sort_values("score_ess", ascending=False, inplace=True)
gene_C_n = Crispr.shape[0]
print("gene_#", gene_C_n, Crispr.shape)
label = np.zeros([gene_C_n])
label[0:top_n] = 1
Crispr['label_ess'] = label
Crispr.sort_values("score_non_ess", ascending=False, inplace=True)
Crispr['label_ness'] = label
print("sum,", sum(Crispr['label_ess'][0:2000]))
print(Crispr)
Crispr = Crispr.loc[genes]

shRNA.sort_values('score_ess', ascending=False, inplace=True)
gene_s_n = shRNA.shape[0]
print("gene#_shRNA", gene_s_n)
label_s = np.zeros([gene_s_n])
label_s[0:top_n] = 1
shRNA['label_ess'] = label_s
shRNA.sort_values('score_non_ess', ascending=False, inplace=True)
shRNA['label_ness'] = label_s
shRNA = shRNA.loc[genes]
print(shRNA)
# print("sum,",sum(shRNA['label_ess'][0:2000]))

f = open("gene_name.txt",'w')
for i in genes:
    f.write(str(i) + '\n')


exp_num = Exp.values
Crispr_num = Crispr.values
print("Crispr_num", Crispr_num)
shRNA_num = shRNA.values
m, n = exp_num.shape

print("m,n", m, n)  # m为基因数，n为cellline数

CRISP_e = np.zeros([n])
CRISP_ne = np.zeros([n])
CRISP_e_S = np.zeros([n])
CRISP_ne_S = np.zeros([n])

shRNA_e = np.zeros([n])
shRNA_ne = np.zeros([n])
shRNA_e_S = np.zeros([n])
shRNA_ne_S = np.zeros([n])

shRNA_essS_exp = []
shRNA_nessS_exp = []
Crispr_essS_exp = []
Crispr_nessS_exp = []
exp_th = 1.3
same_ess_count = 0
shRNA_ess_count = 0
crispr_ess_count = 0

both_ess = []
both_noness = []
only_e_Crispr = []
only_e_shRNA = []
only_ne_Crispr = []
only_ne_shRNA = []
e_Crispr_ne_shRNA = []
ne_Crispr_e_shRNA = []

both_ess_exp = []
both_ness_exp = []
only_e_C_exp = []
only_e_S_exp = []
only_ne_C_exp = []
only_ne_S_exp = []
e_C_ne_S_exp = []
ne_C_e_S_exp = []

marker = np.zeros([m, n*3])
for i in range(n):  # 对于每一个cell line 254
    sort_exp = np.argsort(exp_num[:, i])
    ##########################Calculte the expression of all essential genes in CRISPR and shRNA ######################
    correlation_Crisp_essS_exp = []
    correlation_Crisp_nssS_exp = []
    CRISP_i_sort = np.argsort(Crispr_num[:, i])
    CRISP_ni_sort = np.argsort(-Crispr_num[:, i])
    th_e = Crispr_num[CRISP_i_sort[int(top_rate * m)], i]
    print("the", i, th_e)
    th_ne = Crispr_num[CRISP_ni_sort[int(top_rate * m)], i]
    # 对每一个cellline



    shRNA_i_sort = np.argsort(shRNA_num[:, i])
    th_S_e = shRNA_num[shRNA_i_sort[int(top_rate * m)], i]

    shRNA_i_sort_n = np.argsort(-shRNA_num[:, i])
    th_S_ne = shRNA_num[shRNA_i_sort_n[int(top_rate * m)], i]

    temp_ess_C = []
    temp_ness_C = []
    temp_ess_S = []
    temp_ness_S = []

    for j in range(m):  # 对于每一个gene
        if Crispr_num[j, i] < th_e and Crispr_num[j, n + 2] == 1:  # CRISPR essential
            crispr_ess_count = crispr_ess_count + 1
            # if Crispr_num[j, i] > th_e and Crispr_num[j, n + 3] == 1:
            Crispr_essS_exp.append([Crispr_num[j, i], exp_num[j, i]])
            temp_ess_C.append([j, Crispr_num[j, i], exp_num[j, i]])
            marker[j,3*(i-1)]=1
        elif Crispr_num[j, i] > th_ne and Crispr_num[j, n + 3] == 1:  # CRISPR nonessential
            Crispr_nessS_exp.append([Crispr_num[j, i], exp_num[j, i]])
            temp_ness_C.append([j, Crispr_num[j, i], exp_num[j, i]])
            marker[j,3*(i-1)]=-1
        else:
            marker[j,3*(i-1)]=0


        if shRNA_num[j, i] < th_S_e and shRNA_num[j, n + 2] == 1:  # if shRNA_num[j, i] >th_S_e and shRNA_num[j, n + 3] == 1:  #shRNA essential
            shRNA_ess_count = shRNA_ess_count + 1
            shRNA_essS_exp.append([shRNA_num[j, i], exp_num[j, i]])
            temp_ess_S.append([j, shRNA_num[j, i], exp_num[j, i]])
            marker[j, 3 * (i - 1)+1] = 1
            # non_ess genes in each cell line #shRNA nonessential
        elif shRNA_num[j, i] > th_S_ne and shRNA_num[j, n + 3] == 1:
            shRNA_nessS_exp.append([shRNA_num[j, i], exp_num[j, i]])
            # 里面加入了基因的序列号，shRNA的essential得分，基因表达值
            temp_ness_S.append([j, shRNA_num[j, i], exp_num[j, i]])
            marker[j, 3 * (i - 1)+1] = -1
        else:
            marker[j,3*(i-1)+1]=0

        marker[j, 3 * (i - 1) + 2] = exp_num[j,i]



    # print("shRNA_ess_count", shRNA_ess_count, "Crispr_ess_count", crispr_ess_count)

    temp_ness_S_ID = [i[0] for i in temp_ness_S]
    temp_ess_S_ID = [i[0] for i in temp_ess_S]
    temp_ess_C_ID = [i[0] for i in temp_ess_C]
    temp_ness_C_ID = [i[0] for i in temp_ness_C]

    temp_A = list(set(temp_ess_S_ID).intersection(set(temp_ess_C_ID)))
    both_ess.extend(temp_A)
    both_ess_exp.extend(exp_num[temp_A, i])

    temp_B = list(set(temp_ness_S_ID).intersection(set(temp_ness_C_ID)))
    both_noness.extend(temp_B)
    both_ness_exp.extend(exp_num[temp_B, i])

    temp_G = list(set(temp_ess_C_ID).intersection(set(temp_ness_S_ID)))
    e_Crispr_ne_shRNA.extend(temp_G)
    e_C_ne_S_exp.extend(exp_num[temp_G, i])

    temp_F = list(set(temp_ness_C_ID).intersection(set(temp_ess_S_ID)))
    ne_Crispr_e_shRNA.extend(temp_F)
    ne_C_e_S_exp.extend(exp_num[temp_F, i])

    temp_C = list(set(list(set(temp_ess_C_ID).difference(set(temp_A)))).difference(set(temp_G)))
    only_e_Crispr.extend(temp_C)
    only_e_C_exp.extend(exp_num[temp_C, i])

    temp_D = list(set(list(set(temp_ess_S_ID).difference(set(temp_A)))).difference(set(temp_F)))
    only_e_shRNA.extend(temp_D)
    only_e_S_exp.extend(exp_num[temp_D, i])

    temp_E = list(set(list(set(temp_ness_C_ID).difference(set(temp_B)))).difference(set(temp_F)))
    only_ne_Crispr.extend(temp_E)
    only_ne_C_exp.extend(exp_num[temp_E, i])

    temp_F = list(set(list(set(temp_ness_S_ID).difference(set(temp_B)))).difference(set(temp_G)))
    only_ne_shRNA.extend(temp_F)
    only_ne_S_exp.extend(exp_num[temp_F, i])
    print(i, "shRNA_ess_count", shRNA_ess_count, "Crispr_ess_count", crispr_ess_count)

    # same_ess_count = same_ess_count + len(list(set(temp_shRNA).intersection(set(temp_crispr))))

np.savetxt("marker.csv", marker,delimiter=',')

print("both_ess", len(both_ess), len(list(set(both_ess))), "both_non_ess", len(both_noness),
      len(list(set(both_noness))), "only ess in CRISP", len(only_e_Crispr), len(list(set(only_e_Crispr))),
      "only ess in shRNA", len(only_e_shRNA), len(list(set(only_e_shRNA))),
      "only non ess CRISPR", len(only_ne_Crispr), len(list(set(only_ne_Crispr))), "only_ne_shRNA", len(only_ne_shRNA),
      len(list(set(only_ne_shRNA))),
      "e_Crispr_ne_shRNA", len(e_Crispr_ne_shRNA), len(list(set(e_Crispr_ne_shRNA))), )
print("ne_Crispr_e_shRNA", len(ne_Crispr_e_shRNA), len(list(set(ne_Crispr_e_shRNA))))
print("Draw the ", np.mean(np.array(both_ess_exp)), np.mean(both_ness_exp), np.mean(only_e_C_exp),
      np.mean(only_e_S_exp), np.mean(only_ne_C_exp), np.mean(only_ne_S_exp), np.mean(e_C_ne_S_exp),
      np.mean(ne_C_e_S_exp))

y1 = pd.Series(np.array(both_ess_exp))
y2 = pd.Series(np.array(both_ness_exp))
y3 = pd.Series(np.array(only_e_C_exp))
y4 = pd.Series(np.array(only_e_S_exp))
y5 = pd.Series(np.array(only_ne_C_exp))
y6 = pd.Series(np.array(only_ne_S_exp))
y7 = pd.Series(np.array(e_C_ne_S_exp))
y8 = pd.Series(np.array(ne_C_e_S_exp))
data = pd.DataFrame({"g1": y1, "g2": y3, "g3": y4, "g4": y7, "g5": y8, "g6": y2, "g7": y6, "g8": y5})
# data.boxplot( boxprops = {'color': 'blue'}, whiskerprops={'linestyle':'--','color':'black'}, medianprops = {'color': 'red'}, flierprops = {'marker':'_', 'markerfacecolor':'red'}) # 这里，pandas自己有处理的过程，很方便哦。
# data.boxplot( boxprops = 'b', whiskerprops={'linestyle':'--','color':'black'}, medianprops = 'r', flierprops = {'marker':'_', 'markerfacecolor':'red'}) # 这里，pandas自己有处理的过程，很方便哦。
# flierprops = dict(marker='_', markersize=3, markerfacecolor='r')
############method1#####################
# flierprops = dict(marker='+',  markerfacecolor='green', color = 'red',markersize=3)
# whiskerprops= dict(linestyle='--',linewidth = 1, color='black')
# sns.boxplot(palette="husl", data=data,saturation=1,flierprops=flierprops, whiskerprops=whiskerprops)
############method1#####################

PROPS = {
    'boxprops': {'facecolor': 'none', 'edgecolor': 'blue', 'linewidth': 1},
    'medianprops': {'color': 'red'},
    'whiskerprops': {'linestyle': '--', 'linewidth': 1, 'color': 'black'},
    'capprops': {'color': 'black', 'linewidth': 1},
    'flierprops': {'marker': '+', 'linewidth': 1, 'markeredgecolor': 'r'}
}
sns.boxplot(palette="husl", data=data, saturation=0.8, **PROPS)
plt.ylabel("Gene expression")
plt.xlabel("Groups")  # 我们设置横纵坐标的标题。
plt.show()

##########################Calculte the expression of all essential genes in CRISPR and shRNA ######################
'''
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
    plt.plot(x, y1, lw=1,color = 'teal', label='Crispr_essential')
    plt.plot(x, y11, lw=1, color = 'teal', linestyle='dashed', label='Crispr_non_essential')
    print("y11", y11)
    plt.plot(x, y3, lw=1,color = 'indigo', label='shRNA_essential')
    plt.plot(x, y33, lw=1,  color = 'darkviolet', linestyle='dashed', label='shRNA_non_essential')
    plt.title('top 30% high expression genes')
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
'''
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
