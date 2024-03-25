import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#pd.set_option('display.max_rows', None)

##################################draw the iteration cost value#############################

filename0 = './final/Plato_local/Unsupervised_1000_CRISPR_cost.txt'
Crisper_cost = np.loadtxt(filename0)

filename1 = './final/Plato_local/shRNA_Cost_Unsupervised_1000.txt'
shRNA_cost = np.loadtxt(filename1)

plt.plot(range(Crisper_cost.shape[0]),Crisper_cost[:,0],  lw=1, c='b', label='CRISPR')
# plt.plot(range(Crisper_cost.shape[0]),Crisper_cost[:,1],  lw=1, c='b', label='CRISPR')
plt.plot( range(shRNA_cost.shape[0]),shRNA_cost[:,0], lw=1, c='r', label='shRNA')
# plt.plot( range(shRNA_cost.shape[0]),shRNA_cost[:,1], lw=1, c='r', label='shRNA')

plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
# plt.title('The convergence curves of the objective function ', fontsize = 15)
# plt.xlabel('Iterations',fontsize = 15)
# plt.ylabel('The value of the cost function', fontsize = 15)
plt.title('The convergence curves of matrix U ', fontsize = 15)
plt.xlabel('Iterations',fontsize = 15)
plt.ylabel('The difference of matrix Ut and U(t+1)', fontsize = 15)
plt.legend()
plt.show()



##################################draw the iteration cost value#############################
#
# #the ID of the predicted essential genes
# # filename0 = './shRNA_Essential_genes_top1000_top20per.txt'
# # filename0 = './final/ess_CRISPR_top1500.txt'
# filename0 = './final/shRNA/ess_shRNA_top1000.txt'
#
# #the # of essential genes in each celllines by 10%
# # filename1 = './shRNA_essential.txt'
# filename1 = './final/shRNA/shRNA_essential_top0.1_ess_cellline_sum.txt'
#
# #the # of non-Essential genes in each cellines by 10%
# # filename1 = './shRNA_essential.txt'
# filename2 = './final/shRNA/shRNA_non_essential_top0.1_noness_cellline_sum.txt'
#
#
# # filename2 = './shRNA_non_essential.txt'
# id = pd.read_table(filename0,sep='\t',header=None)
# esse = pd.read_table(filename1,sep='\t',header=None)
# n_esse = pd.read_table(filename2,sep='\t',header=None)
#
# df = pd.concat([id, esse, n_esse], axis=1)
# df.columns = ['id', 'essential', 'non_essential']
# df.sort_values(['id', 'non_essential'], ascending=False, inplace= True)
#
#
# y1 = df['non_essential'].values.tolist()
# y2 = df['essential'].values.tolist()
# y0 = (df['id'].values*2).tolist()
#
# print("ID",sum(y0), sum(y2[0:1468])/1468 )
# def plot_line(df, y1, y2, y0):
#     x = range(df.shape[0])
#     plt.plot(x, y2, lw=1, c='#d62728', label='total cell lines with essential gene')  #
#     plt.plot(x, y1, lw=1, c='b', label='total cell lines with non_essential gene')  ##1f77b4
#
#     x1=[]
#     y11=[]
#     x2=[]
#     y12=[]
#     for i in range(df.shape[0]):
#         if y0[i]>0:
#             x1.append(i)
#             y11.append(y0[i])
#         else:
#             x2.append(i)
#             y12.append(y0[i])
#     plt.plot(x1, y11, lw=1, c='b', label='Essential genes')
#     plt.plot(x2, y12, lw=1, c='r', label='other genes')
#     # plt.scatter(x, y0, c='#2ca02c', s= 0.5, label='0_1')
#
#
#     plt.title('(C) Essential common genes from shRNA ')
#     plt.xlabel('Gene')
#     plt.ylabel('Cell line number')
#     #plt.xticks(x, genes, rotation=45)
#     plt.legend()
#     plt.show()
#
#
#
#
#
# plot_line(df, y1, y2, y0)
