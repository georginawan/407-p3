import numpy as np
import xlrd
import scipy.stats as st
# Initialize effect_size, control_mean, control_sd
effect_size, sample_size, control_mean, control_sd = 0.05, 50, 1, 0.5

ge_file = open("/Users/yuwan/Downloads/TCGA-LUAD&LUSC/luad_data_new.txt")
ge_data = ge_file.readlines()
Case_id = ge_data[0][:-2].split(",")
ge_mat = []
for i in range(2,len(ge_data)):
    line = [float(i) for i in ge_data[i][:-2].split(",")[1:]]
    ge_mat.append(line)
ge_mat = np.array(ge_mat)


clinical_file = xlrd.open_workbook(r"/Users/yuwan/Downloads/luad_clinical_cbio.xls")
sheet1 = clinical_file.sheet_by_index(0)
feature = sheet1.row_values(0,start_colx=0,end_colx=None)
ddt_index = feature.index('Overall Survival (Months)')
status_index = feature.index('Overall Survival Status')
Death_days_to = sheet1.col_values(ddt_index,start_rowx=1,end_rowx=None)
Survival_status = sheet1.col_values(status_index,start_rowx=1,end_rowx=None)
samples_cbio = sheet1.col_values(2,start_rowx=1, end_rowx=None)
year = 2
label_list = []
for i in range(len(Case_id)):        
    sample_index_cbio = samples_cbio.index(Case_id[i])
    label = float(Death_days_to[sample_index_cbio])
    label = round(label)
    survival_status = (Survival_status[sample_index_cbio])
    #sample_index = samples.index(Case_id[i])
    #contact_days = float(last_contact_days_to[sample_index])
    if survival_status == "0:LIVING":
        if label >= 12*year:
            label_list.append(1)
        else:
            pass

    else:
        if label >= 12*year:
            label_list.append(1)            
        else:
            label_list.append(2)
    
state1,state2 = [],[]
for i in range(len(label_list)):
    if label_list[i] == 1:
        state1.append(ge_mat[:,i])
    else:
        state2.append(ge_mat[:,i])
state1 = np.array(state1)
state2 = np.array(state2)
n_1,n_gene = np.shape(state1)
n_2,n_gene = np.shape(state2)


def power_analysis(sample_size = 10,sims = 10,data1=[],data2=[]):
    data1_ = data1[:sample_size,:]
    data2_ = data2[:sample_size,:]
    power_list = []
    # Keep incrementing sample size by 10 till we reach required power
    while True:
        p_count = 0
        for i in range(sims):
            mean1 = np.mean(data1_[:,i])
            std1 = np.std(data1_[:,i])
            mean2 = np.mean(data2_[:,i])
            std2 = np.std(data2_[:,i])
            control_time_spent = np.random.normal(loc=mean1, scale=std1, size=(sample_size,1))
            treatment_time_spent = np.random.normal(loc=mean2, scale=std2, size=(sample_size,1))
            t, p = st.ttest_ind(treatment_time_spent, control_time_spent)
            if p < 0.05:
                p_count += 1
    
        print(sample_size,p_count/sims)
        power_list.append(p_count/sims)
        if p_count/sims >= 0.8: 
            print(power_list)
            break
        elif sample_size>=1000:
            print(sample_size,p_count/sims)
            print(power_list)
            break
        else:
            sample_size += 50

    print("For 80% power, sample size required = {}".format(sample_size))
    
power_analysis(10,n_gene,state1,state2)
ge_file.close()
