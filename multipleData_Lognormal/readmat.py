import scipy.io

#results_Wass_50_new_2.mat: [1,10,35,65]=>[10, 20, 25, 30]
#results_Wass_500_new_2.mat: [1,10,35,65]=>[10, 20, 25, 30]
#results_Wd_50_new_2.mat: [1,10,35,65]=>[10, 20, 42, 60]
#results_Wd_500_new_2.mat: [1,10,35,65]=>[10, 20, 43, 60]

# filename='results_Wass_50_new_2.mat'
# filename='results_Wass_500_new_2.mat'
# filename='results_Wd_50_new_2.mat'
filename='results_Wd_500_new_2.mat'


mat = scipy.io.loadmat('./sanjula_results/'+filename)
#results=mat['Wass_DRO_zobj']
results=mat['Wd_DRO_zobj']

k=60
print(max(results[:,k-1:k][0:30]))
print(min(results[:,k-1:k][0:30]))