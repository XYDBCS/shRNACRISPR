#From the predicted score to get the top 1000 essential genes ID and non_essential gene ID, besides get the intersect essentail genes in each cell lines

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# A =np.array([[1,2,3], [4,2,1]])
# B = np.array([[2,1,3], [2,4,1]])
# print("A", "B", A, B)
# print(np.linalg.norm(A-B))

#Essential score
# filename1 = "./data/Essentiality shRNA 254 cell lines.xlsx"
filename1 = "./data/Essentiality CRISPR 254 cell lines.xlsx"
Crispr = pd.read_excel(filename1, index_col=0) #Column (0-indexed) to use as the row labels of the DataFrame.
print("shRNA shape", Crispr.shape)
genes1 = Crispr.index.tolist()  #get genes name

filename0 = './data/Gene Expression 530 Cell lines.xlsx'
Exp = pd.read_excel(filename0, index_col=0)
genes_exp = Exp.index.tolist()  # genes with expression info



# #Essential gene #the known essential genes
# # filename0 = "./data/shRNA_training_essentials.xlsx"
# filename0 = "./data/Essentials_CRISPR.xlsx"
# Ess = pd.read_excel(filename0, index_col=0) #the first one in each row as the index
# # Ess.sort_values('Gene ID', inplace=True) #sort by gene ID and all the other elements follow it
# print("ESS.shape",Ess.shape)
# genes0 = Ess.index.tolist()   #genes name
#
# #add one more colum for mark the essential genes
# Crispr['essential'] = Crispr.index.isin(Ess.index)
# print("sum_essential",sum(Crispr.index.isin(Ess.index)))
# print(Ess.loc[~Ess.index.isin(Crispr.index)]) #Print the one that not in CRISPR score

#choose how many persent as the essential or non-essential genes
threshold = 0.1

#Change the essential score to essential information with -1(non-essential), 0(irrelevant), 1(essential)
# print("dataFram_Crisp", Crispr)
# print("Crispr",Crispr.values)
ess_scor_n = Crispr.values
#
# #get the standard essential genes
# ess_ID = ess_scor_n[:,-1].astype(int)#one row
#
# print("ess_ID",ess_ID)

print("ess_score_n", ess_scor_n)
#delete the first coulme gene ID and the essential infor, keep the score information with one row is one gene, one colum is a cellline
# ess_scor_n = np.delete(ess_scor_n,-1, axis = 1)
ess_scor_n = np.delete(ess_scor_n, [0], axis = 1)#delete the first column of gene ID


#m is the number of gene, n is the # of cell line
m,n = ess_scor_n.shape


ess_int_n = np.zeros([m,n])
#along x-axis, sort each colum
sort_ess_index = np.argsort(ess_scor_n, axis = 0)
#print("sort", sort_ess_index[:,0], max(sort_ess_index[:,0]))
for i in range (n):
    ess_th = m*threshold
    non_ess_th = m*(1-threshold)
    for j in range(m):
        #top small % as essential gene with 1
        if sort_ess_index[j,i] < ess_th:
        # if ess_scor_n[j, i] < 0:
            ess_int_n[j,i] = 1
        #top last % large value as non-essential gene with -1
        elif sort_ess_index[j,i] >non_ess_th:
        # elif ess_scor_n[j, i] > 0:
            ess_int_n[j,i] = -1

essential = np.zeros(m)
non_essential = np.zeros(m)
for i in range(m):
    essential[i] = np.sum(ess_int_n[i,:]>=1) #根据top m*threshold 的基因为essential，每个基因在多少cell line里是essential的
    non_essential[i] = np.sum(ess_int_n[i,:]<=-1)#根据降序top m*threshold 的基因为non_essential，每个基因在多少cell line里是non_essential的

print("CRISPR_Score_essential, CRISPR_Score_non_essential",essential, non_essential)

# np.savetxt("./final/shRNA_essential_top0.1_ess_cellline_sum.txt", essential)
# np.savetxt("./final/shRNA_non_essential_top0.1_noness_cellline_sum.txt", non_essential)
# np.savetxt("shRNA_ess_ID_by0.txt", ess_ID)

# filename2 = "./Plato_local/supervised_prediction score_0.15_0.0001_norm_256_6000.txt"
# filename2 = "./supervised_prediction score_RQ_0.15_0.0001_6000_256.txt"
# filename2 = "./shRNA_Semi_Del_prediction score_0.1_0.0001_norm_256_6000.txt"
filename2 = "./final/Plato_local/Unsupervised_1000_CRISPR.txt"
# filename2 = "./final/Plato_local/shRNA_Unsupervised_1000.txt"
predict = np.loadtxt(filename2)

Predict_rank = np.argsort(-predict, axis = 0)#rank by descenting, big-small
top_n = 1000
ess_overlap = 0
# print("rank score",predict[Predict_rank[1900:2000,0]])
k = 0#when k = 0, its essential


# filename3 = "./Plato_local/unsupervised_prediction score_0.15_0.0001_norm_6000_256.txt"
# filename3 = "./Unsupervised_prediction score_RQ_0.15_0.0001_256_6000.txt"
# filename3 = "./shRNA_unsupervised_Del_prediction score_0.1_0.0001_norm_256_6000.txt"
#
# predict2 = np.loadtxt(filename3)
# Predict_rank2 = np.argsort(-predict2, axis = 0)
# intersection = np.intersect1d(Predict_rank[0:top_n,k],Predict_rank2[0:top_n,k])
# print("intersection",intersection, intersection.shape)
intersection = Predict_rank[0:top_n,k]#the ID of top n essential genes from unsupervised learning
ess_shRNA_score = np.zeros([top_n])
ess_shRNA_score = predict[intersection,0]
print("essential", intersection[0:10], predict[intersection[0],0], predict[8245,0],ess_shRNA_score[0:3])
# np.savetxt("ess_shRNA_score.txt",ess_shRNA_score)
np.savetxt("ess_CRISPR_score.txt",ess_shRNA_score)




non_ess_rank = Predict_rank[0:top_n, 2]#The ID of the top n non_essential genes
noness_shRNA_score = np.zeros([top_n])
noness_shRNA_score = predict[non_ess_rank,2]
print("non_ess", non_ess_rank[0:10], predict[non_ess_rank[0],2],noness_shRNA_score[0:3])
# np.savetxt("noness_shRNA_score.txt", noness_shRNA_score)
np.savetxt("noness_CRISPR_score.txt", noness_shRNA_score)

ess = np.zeros(m)
ess_ID = intersection

non_ess = np.zeros(m)
noness_ID = non_ess_rank

shRNA_predicted_1000_common_essential = []
shRNA_predicted_1000_common_nonessential = []
for i in range(intersection.shape[0]):
    ess[intersection[i]]=1
    shRNA_predicted_1000_common_essential.append(genes1[intersection[i]])
    non_ess[non_ess_rank[i]] = 1
    shRNA_predicted_1000_common_nonessential.append(genes1[non_ess_rank[i]])


'''

# np.savetxt("shRNA_non_Essential_genes_ID_top1500_top20per.txt", intersection)
# np.savetxt("shRNA_non_Essential_genes_top1500_top20per.txt", ess)
files = open('./final/shRNA/ess_shRNA_top1000_ID_new.txt', 'w')
files.write(str(shRNA_predicted_1000_common_essential))
files.close()

filesn = open('./final/shRNA/noness_shRNA_top1000_ID_new.txt', 'w')
filesn.write(str(shRNA_predicted_1000_common_nonessential))
filesn.close()

print("ess", ess.shape) #ess 里面只有0,1， 当为1时表示一个ess gene， 0表示一个non_ess gene
# np.savetxt("./final/shRNA/ess_shRNA_top1000.txt", ess)
# np.savetxt("./final/shRNA/ess_shRNA_top1000_ID.txt", ess_ID)
# np.savetxt("./final/shRNA/noness_shRNA_top1000.txt", non_ess)
# np.savetxt("./final/shRNA/noness_shRNA_top1000_ID.txt", noness_ID)

filename1 = "./data/Essentiality CRISPR 254 cell lines.xlsx"
Crispr = pd.read_excel(filename1, index_col=0)
cell_lines = Crispr.columns.tolist()  # cell line name from CRISPR
Exp = Exp.reindex(cell_lines, axis='columns')


Exp1=Exp
genes_1000 = list(set(genes_exp).intersection(set(shRNA_predicted_1000_common_essential)))
genes_non_1000 =  list(set(genes_exp).intersection(set(shRNA_predicted_1000_common_nonessential)))
Exp_ess_CRISPR = Exp.loc[genes_1000].drop(columns=['GeneID'])
Exp_ess_CRISPR = np.array(Exp_ess_CRISPR)
Exp_noness_CRISPR = Exp1.loc[genes_non_1000].drop(columns=['GeneID'])
Exp_noness_CRISPR = np.array(Exp_noness_CRISPR)


print(Exp_noness_CRISPR.shape, Exp_ess_CRISPR.shape)


np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=np.inf)

file = open('./final/shRNA/predicted_ess_shRNA_1000.txt', 'w')
file.write(str(Exp_ess_CRISPR))
file.close()

# np.set_printoptions(threshold=np.inf)
file=open('./final/shRNA/predicted_noness_shRNA_1000.txt','w')
file.write(str(Exp_noness_CRISPR))
file.close()


# for i in range(top_n):
#     if ess_ID[Predict_rank[i,k]] ==1:
#         ess_overlap = ess_overlap+1
# print("overlap with essential gene", ess_overlap)



'''