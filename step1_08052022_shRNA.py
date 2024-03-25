import pandas as pd
#import matplotlib.pyplot as plt
import numpy as np

#Essential score
filename1 = "./data/Essentiality shRNA 254 cell lines.xlsx"
Crispr = pd.read_excel(filename1, index_col=0,engine='openpyxl') #Column (0-indexed) to use as the row labels of the DataFrame.
#Crispr = pd.read_excel(filename1, index_col=0) #Column (0-indexed) to use as the row labels of the DataFrame.
print("shRNA shape", Crispr.shape)
genes1 = Crispr.index.tolist()  #get genes name

#Essential gene
filename0 = "./data/shRNA_training_essentials.xlsx"
Ess =  pd.read_excel(filename0, index_col=0,engine='openpyxl')  #the first one in each row as the index
#Ess = pd.read_excel(filename0, index_col=0) #the first one in each row as the index
# Ess.sort_values('Gene ID', inplace=True) #sort by gene ID and all the other elements follow it
print("ESS.shape",Ess)
genes0 = Ess.index.tolist()   #genes name
print(genes0)

#non_essential gene
filename2 = "./data/shRNA_training_nonessential.xlsx"
Ess2 =  pd.read_excel(filename2, index_col=0,engine='openpyxl')  #the first one in each row as the index
#Ess = pd.read_excel(filename0, index_col=0) #the first one in each row as the index
# Ess.sort_values('Gene ID', inplace=True) #sort by gene ID and all the other elements follow it
print("ESS2.shape",Ess2)
genes2 = Ess2.index.tolist()   #genes name
print(genes2)


#add one more colum for mark the essential genes
Crispr['essential'] = Crispr.index.isin(Ess.index)
Crispr['non_essential'] = Crispr.index.isin(Ess2.index)
print("sum_essential",sum(Crispr.index.isin(Ess.index)), sum(Crispr.index.isin(Ess2.index)))
# print(Ess.loc[~Ess.index.isin(Crispr.index)]) #Print the one that not in CRISPR score

#choose how many persent as the essential or non-essential genes
threshold = 0.1

#Change the essential score to essential information with -1(non-essential), 0(irrelevant), 1(essential)
print("dataFram_Crisp", Crispr)
print("Crispr",Crispr.values)
ess_scor_n = Crispr.values

#get the standard essential genes
ess_ID = ess_scor_n[:,-2].astype(int)#one row
print("ess_ID", sum(ess_ID))

non_ess_ID = ess_scor_n[:,-1].astype(int)
print("non_ess", sum(non_ess_ID))

#delete the first coulme gene ID and the essential infor, keep the score information with one row is one gene, one colum is a cellline
ess_scor_n = np.delete(ess_scor_n,-1, axis = 1)
ess_scor_n = np.delete(ess_scor_n,-1, axis = 1)
ess_scor_n = np.delete(ess_scor_n, [0], axis = 1)

#m is the number of gene, n is the # of cell line
m,n = ess_scor_n.shape
ess_int_n = np.zeros([m,n])

#set an matrix with essential groups (row) and gene (colum)
ess_gene = np.zeros([3*n, m])

#along x-axis, sort each colum
sort_ess_index = np.argsort(ess_scor_n, axis = 0)
#print("sort", sort_ess_index[:,0], max(sort_ess_index[:,0]))
for i in range (n):
    ess_th = m*threshold
    non_ess_th = m*(1-threshold)
    for j in range(m):
        #top small % as essential gene with 1
        if sort_ess_index[j,i] < ess_th:
            ess_int_n[j,i] = 1
            ess_gene[i*3+0, j] = 1
        #top last % large value as non-essential gene with -1
        elif sort_ess_index[j,i] >non_ess_th:
            ess_int_n[j,i] = -1
            ess_gene[i*3+2, j] =1
        else:
            ess_gene[i*3+1, j] =1




#calculate the cellLine # with essential genes and non_essential genes
essential = np.zeros(m)
non_essential = np.zeros(m)
for i in range(m):
    essential[i] = np.sum(ess_int_n[i,:]>=1)
    non_essential[i] = np.sum(ess_int_n[i,:]<=-1)

#only use the top 20% essential genes
real_essential = np.zeros(m)
real_nonessential = np.zeros(m)
real_ess_id = np.zeros(m)
real_nonessential_id = np.zeros(m)

for i in range(m):
    if ess_ID[i]==1:
        real_essential[i] = essential[i]
    if non_ess_ID[i] ==1:
        real_nonessential[i] = non_essential[i]


top_n = int(sum(ess_ID)*0.2)
top_n_noness = int(sum(non_ess_ID)*0.2)

sort_real_ess_index = np.argsort(-real_essential)
for i in range(top_n):
    real_ess_id[sort_real_ess_index[i]] = 1
print("top 20%", real_essential[sort_real_ess_index])
print("Average_essential", sum(real_essential[:])/sum(ess_ID[:]))
ess_ID = real_ess_id

sort_real_noness_index = np.argsort(-real_nonessential)
for i in range(top_n_noness):
    real_nonessential_id [sort_real_noness_index[i]] = 1
print("top 20%_nonessential", real_nonessential[sort_real_noness_index])
print("Average_nonessential", sum(real_nonessential[:])/sum(non_ess_ID[:]))
non_ess_ID = real_nonessential_id
print("after_non_essential_ID", sum(non_ess_ID))





#print("sum", sum(ess_gene[0,:]))
Dn_tep = np.zeros([m])
Dn_tep[:] = n
Dn = np.diag(Dn_tep)
#print("Dn", Dn)

Dv_diag_array = np.zeros([3*n])
for i in range(n):
    Dv_diag_array[i*3-3] = ess_th
    Dv_diag_array[i*3-2] = m-2*ess_th
    Dv_diag_array[i*3-1] = ess_th
Dv = np.diag(Dv_diag_array)
#print("dv", Dv)

#F is the matrix with essential genes
F = np.zeros([m,3])
F[:,0] = ess_ID.T
F[:,2] = non_ess_ID.T
#print(F)

Y = np.zeros([n*3, 3])
for i in range(n):
    Y[i*3-3, 0] = 1
    Y[i*3-2, 1] = 1
    Y[i*3-1, 2] = 1
# print("Y", Y)

Kv = np.eye(n*3, k = 0)
# print("Kv", Kv)
diag_ID = ess_ID+non_ess_ID
Hn = np.diag(diag_ID)
# print("Hn", Hn)

#initial the Ut by random value
# Ut1 = np.zeros([m,3])
# for i in range(m):
#     Ut1[i,:] = np.random.dirichlet(np.ones(3),  size=1)
# # print("Ut1", Ut1)

a = 6000
b = 256
gene_n = m
cell_v = n*3
diff = 100
diff_th = 0.0001

def sum1_norm(a):
    m,n = a.shape
    b = np.zeros([m,n])
    for i in range(m):
        b[i,:] = a [i,:]/sum(a[i,:])
    return b

# Qt = np.linalg.inv(Dv + a * Kv).dot((ess_gene).dot(Ut1) + a * Kv.dot(Y))
# Qt = sum1_norm(Qt)
Qt = Y
Ut1 = np.linalg.inv(Dn).dot(ess_gene.T).dot(Qt)
print("ut1", Ut1)
I = 0
while diff > diff_th and I<100:
    # Qt = np.linalg.inv(Dv + a * Kv).dot((ess_gene).dot(Ut1) + a * Kv.dot(Y))
    # print("Qt", Qt)
    Qt = np.linalg.inv(Dv + a * Kv).dot((ess_gene).dot(Ut1) + a * Kv.dot(Y))
    # Qt = sum1_norm(Qt)
    print("Qt", Qt)
    #Ut = np.linalg.inv(Dn + b * Hn).dot((ess_gene.T).dot(Qt) + b * Hn.dot(F))
    Ut = np.linalg.inv(Dn).dot(ess_gene.T).dot(Qt)
    # Ut = sum1_norm(Ut)
    print("Ut", Ut)
    diff = np.linalg.norm(Ut-Ut1)
    Ut1 = Ut
    I = I + 1
    print("i", I, diff)

np.savetxt("shRNA_unsupervised_Del_prediction score_0.1_0.0001_norm_256_6000.txt", Ut)
