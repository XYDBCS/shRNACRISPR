import pandas as pd
import matplotlib.pyplot as plt
# import seaborn as sns
import numpy as np
import seaborn as sns

#######################################task 1: draw the distributation of genes and task 2: distributation of gene expression ###################################
####################################################################################################
top_rate = 0.1
top_rate_pre = 0.1
top_n = int(top_rate_pre*15029)
top_exp_rate = top_rate
filename1 = "./data/Essentiality CRISPR 254 cell lines.xlsx"  #254个 cell line在CRIPSR基因敲除后的值
# filename11 = "./Unsupervised_prediction score_RQ_0.1_0.0001_256_6000.txt"
filename11 ="./unsupervised_prediction score_0.15_0.0001_norm_6000_256.txt"  #预测的得分,没有带geneID
# filename11 ="./final/Plato_local/Unsupervised_1000_CRISPR.txt"
# Crispr = pd.read_excel(filename1, index_col=0,engine='openpyxl')
Crispr = pd.read_excel(filename1, index_col=0)#以第一列作为行索引读数据
cell_lines = Crispr.columns.tolist()#得到CRispr的cell line
# print("cell_lines",cell_lines)

pred_Crispr = np.loadtxt(filename11)
Crispr = Crispr.reindex(cell_lines, axis='columns')


Crispr['score_ess'] = pred_Crispr[:,0]
Crispr['score_non_ess'] = pred_Crispr[:,2]
gene_crispr = Crispr.index.tolist()


filename2 = "./data/Essentiality shRNA 254 cell lines.xlsx"
filename22 = "./shRNA_unsupervised_Del_prediction score_0.1_0.0001_norm_256_6000.txt"
# filename22 ="./final/Plato_local/shRNA_Unsupervised_1000.txt"

# shRNA = pd.read_excel(filename2, index_col=0,engine='openpyxl')
shRNA = pd.read_excel(filename2, index_col=0)
shRNA = shRNA.fillna(0)
pred_shRNA = np.loadtxt(filename22)
shRNA = shRNA.reindex(cell_lines, axis='columns')
shRNA['score_ess'] = pred_shRNA[:,0]
shRNA['score_non_ess'] = pred_shRNA[:,2]


#获取基因表达数据
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
print("m,n", m, n)
CRISP_e = np.zeros([n])
CRISP_ne = np.zeros([n])
CRISP_e_S = np.zeros([n])
CRISP_ne_S = np.zeros([n])

shRNA_e = np.zeros([n])
shRNA_ne = np.zeros([n])
shRNA_e_S = np.zeros([n])
shRNA_ne_S = np.zeros([n])

corr_non_shRNA = []
corr_shRNA = []
corre_Crispr = []
corrn_Crispr = []

shRNA_essS_exp = []
Crispr_essS_exp = []
exp_th=1.8
same_ess_count =0
shRNA_ess_count =0
crispr_ess_count = 0

# n is cell line number
Crispr_low_exp_ess = []
shRNA_low_exp_ess = []

Crispr_all_ess = []
shRNA_all_ess = []
for i in range(n):
    sort_exp = np.argsort(exp_num[:,i])
 ##########################Calculte the expression of all essential genes in CRISPR and shRNA ######################
    correlation_Crisp_essS_exp = []
    correlation_Crisp_nssS_exp = []
    CRISP_i_sort = np.argsort(Crispr_num[:, i])
    # CRISP_i_sort = np.argsort(-Crispr_num[:, i])
    th_e = Crispr_num[CRISP_i_sort[int(top_rate * m)], i]
    # 对每一个cellline
    both_ess = []
    both_noness = []
    only_e_Crispr = []
    only_e_shRNA = []
    only_ne_Crispr =[]
    only_ne_shRNA = []
    e_Crispr_ne_shRNA = []
    ne_Crispr_e_shRNA = []

    temp_crispr = []
    # m is the number of genes
    for j in range(m):
        if  Crispr_num[j, i] < th_e and Crispr_num[j, n+2]==1 and Crispr_num[j, i] > -3.5:

        # if Crispr_num[j, i] > th_e and Crispr_num[j, n + 3] == 1 and Crispr_num[j, i]<2:
        # if Crispr_num[j, i] > th_e and Crispr_num[j, n + 3] == 1:
            Crispr_essS_exp.append([Crispr_num[j, i], exp_num[j, i]])
            Crispr_all_ess.append(j)
            if exp_num[j,i]<exp_th:
                temp_crispr.append(j)
                Crispr_low_exp_ess.append(j)
                crispr_ess_count = crispr_ess_count + 1


    temp_shRNA =[]
    # shRNA_i_sort = np.argsort(-shRNA_num[:, i])
    shRNA_i_sort = np.argsort(shRNA_num[:, i])
    th_S_e = shRNA_num[shRNA_i_sort[int(top_rate * m)], i]
    for j in range(m):
         if shRNA_num[j,i] < th_S_e and shRNA_num[j, n+2] ==1 and shRNA_num[j, i] > -3.5:

        # if shRNA_num[j, i] >th_S_e and shRNA_num[j, n + 3] == 1 and shRNA_num[j, i]<2:
            shRNA_essS_exp.append([shRNA_num[j,i], exp_num[j,i]])
            shRNA_all_ess.append(j)
            if exp_num[j,i]< exp_th:
                temp_shRNA.append(j)
                shRNA_low_exp_ess.append(j)
                shRNA_ess_count = shRNA_ess_count + 1


    same_ess_count = same_ess_count + len(list(set(temp_shRNA).intersection(set(temp_crispr))))# count how many essential genes are identified by both method

Unique_CRISPR = 0
for i in Crispr_low_exp_ess:
    Unique_CRISPR = Unique_CRISPR + 1/Crispr_low_exp_ess.count(i)

Unique_shRNA = 0
for i in shRNA_low_exp_ess:
    Unique_shRNA = Unique_shRNA + 1/shRNA_low_exp_ess.count(i)

Unique_CRISPR_all = 0
for i in Crispr_all_ess:
    Unique_CRISPR_all = Unique_CRISPR_all + 1/Crispr_all_ess.count(i)

Unique_shRNA_all = 0
for i in shRNA_all_ess:
    Unique_shRNA_all = Unique_shRNA_all + 1/shRNA_all_ess.count(i)


print("the common genes in both shRNA and CRISPR", len(list(set(Crispr_all_ess).intersection(set(shRNA_all_ess)))))
print("CRISPR",len(list(set(Crispr_all_ess))))
print("shRNA", len(list(set(shRNA_all_ess))))
Crispr_essS_exp = np.array(Crispr_essS_exp)
shRNA_essS_exp = np.array(shRNA_essS_exp)
print("Crispr_list", Crispr_essS_exp.shape, Crispr_essS_exp)
print("shRNA_list", shRNA_essS_exp.shape, shRNA_essS_exp)


df_shRNA = pd.DataFrame(shRNA_essS_exp,index=None,columns = ['score','expression'])
df_Crispr = pd.DataFrame(Crispr_essS_exp, index=None, columns=['score', 'expression'])
df_shRNA['label'] = 'shRNA:44563 (1322 unique genes)'
df_Crispr['label'] = 'CRISPR:41167 (506 unique genes)'
# df_shRNA['label'] = 'shRNA:44563'
# df_Crispr['label'] = 'Crispr:41167'
# df_shRNA['label'] = 'shRNA:44996'
# df_Crispr['label'] = 'Crispr:36698'

print("low expression essential genes", np.sum(shRNA_essS_exp[:,1]<exp_th), np.sum(Crispr_essS_exp[:,1]<exp_th), Unique_shRNA_all, Unique_CRISPR_all )
print("same ess in low expression", same_ess_count, shRNA_ess_count, crispr_ess_count, Unique_shRNA, Unique_CRISPR)
data = pd.concat([df_shRNA, df_Crispr],sort=False,axis=0,ignore_index=True)
data = data.loc[~(data==0).all(axis=1)]

pal=sns.color_palette("CMRmap")

sns.set(font_scale=1.2)#将字体放大指定的倍数
h =sns.jointplot(x='score',y='expression',data=data,hue='label', s = 4,
                 # joint_kws= dict(alpha = 0.6, fill=True),
                 joint_kws= dict(alpha = 0.6, fill=True),
                 marginal_kws=dict(fill = True),
                 # kind="kde", cmap = pal,
                 kind="scatter", cmap = pal,
                 # palette = pal[1:3]
                 palette={
                     # 'shRNA:44996 (1415 unique genes)': 'blueviolet',
                     # 'CRISPR:36698 (1274 unique genes)': 'firebrick'
                     'shRNA:44563 (1322 unique genes)': 'blueviolet',
                     'CRISPR:41167 (506 unique genes)': 'firebrick'
                    })
    # 'Crispr': 'firebrick'
    # 'shRNA': 'blueviolet',
    # 'shRNA:' + str(shRNA_essS_exp.shape[0])+' (1415 unique genes)': 'blueviolet',
    # 'Crispr:' + str(Crispr_essS_exp.shape[0])+' (1274 unique genes)': 'firebrick'

    # 'shRNA:'+str(shRNA_essS_exp.shape[0]): 'crimson',
    # 'Crispr:'+str(Crispr_essS_exp.shape[0]): 'blue'}

# sns.move_legend(h, "lower left")
plt.legend("lower left")
# h.set_axis_labels('Essential score', 'Gene expression', fontsize=18)

plt.savefig("essential-scatter2.png")
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