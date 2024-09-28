import scipy.io
import numpy as np
from sklearn.mixture import GaussianMixture
import subprocess
import concurrent.futures
import pickle

mat = scipy.io.loadmat('samples_50.mat')
#mat = scipy.io.loadmat('samples_500.mat')
print(len(mat['xi']))

def run_gams(gams_file):
    print(f"Task starts: {gams_file}")  # 打印任务开始信息
    result = subprocess.run(['/mnt/d/gams/47/gams.exe', gams_file], capture_output=True, text=True)
    print(f"Task ends: {gams_file}")  # 打印任务结束信息
    return result.stdout

resultOBJ={}

#for sampleOrder in range(0,len(mat['xi'])):
for sampleOrder in range(0,50):
    resultOBJ[sampleOrder]={}
    sample=np.hstack(mat['xi'][sampleOrder])
    gm = GaussianMixture(n_components=3, covariance_type="diag", random_state=0).fit(sample)
    # print(gm.weights_)
    # print(gm.means_)
    # print(gm.covariances_)
    print('########################',sampleOrder)

    with open('wUH.txt', 'w') as f:
        for i, value in enumerate(gm.weights_, 1):
            f.write(f"wUH('UH{i}') = {value:.3f};\n")

    with open('wLH.txt', 'w') as f:
        for i, value in enumerate(gm.weights_, 1):
            f.write(f"wLH('LH{i}') = {value:.3f};\n")

    # Open the file to write the data
    with open('miuUH.txt', 'w') as f:
        for i, row in enumerate(gm.means_, 1):
            for j, value in enumerate(row, 1):
                f.write(f"miuUH('UH{i}','{j}')={value:.6f};\n")

    # Open the file to write the data
    with open('alphaUH.txt', 'w') as f:
        for i, row in enumerate(gm.covariances_, 1):
            for j, value in enumerate(row, 1):
                f.write(f"alphaUH('UH{i}','{j}')={value:.6f};\n")
    
    gams_files = ['LULinearPortfolioShortLoopMax_1.gms', 'LULinearPortfolioShortLoopMin_1.gms',
                  'LULinearPortfolioShortLoopMax_10.gms', 'LULinearPortfolioShortLoopMin_10.gms',
                  'LULinearPortfolioShortLoopMax_35.gms', 'LULinearPortfolioShortLoopMin_35.gms',
                  'LULinearPortfolioShortLoopMax_65.gms', 'LULinearPortfolioShortLoopMin_65.gms']

    with concurrent.futures.ProcessPoolExecutor() as executor:
        # 提交多个任务并行运行
        results = list(executor.map(run_gams, gams_files))

    file_path = './resultsMAX_1.txt'
    # Open and read the file
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # Process the lines to extract the numerical values after the '=' sign
    resultsMAX_1 = []
    for line in lines:
        if '=' in line:
            # Split the line by '=' and strip leading/trailing spaces
            value = line.split('=')[1].strip()
            # Convert the value to a float and append to the resultsMAX_1 list
            resultsMAX_1.append(float(value))

    file_path = './resultsMIN_1.txt'
    # Open and read the file
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # Process the lines to extract the numerical values after the '=' sign
    resultsMIN_1 = []
    for line in lines:
        if '=' in line:
            # Split the line by '=' and strip leading/trailing spaces
            value = line.split('=')[1].strip()
            # Convert the value to a float and append to the resultsMIN_1 list
            resultsMIN_1.append(float(value))

    resultOBJ[sampleOrder]['1_upper']=resultsMAX_1
    resultOBJ[sampleOrder]['1_lower']=resultsMIN_1

######################## 10 ##############################
    file_path = './resultsMAX_10.txt'
    # Open and read the file
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # Process the lines to extract the numerical values after the '=' sign
    resultsMAX_10 = []
    for line in lines:
        if '=' in line:
            # Split the line by '=' and strip leading/trailing spaces
            value = line.split('=')[1].strip()
            # Convert the value to a float and append to the resultsMAX_10 list
            resultsMAX_10.append(float(value))

    file_path = './resultsMIN_10.txt'
    # Open and read the file
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # Process the lines to extract the numerical values after the '=' sign
    resultsMIN_10 = []
    for line in lines:
        if '=' in line:
            # Split the line by '=' and strip leading/trailing spaces
            value = line.split('=')[1].strip()
            # Convert the value to a float and append to the resultsMIN_10 list
            resultsMIN_10.append(float(value))

    resultOBJ[sampleOrder]['10_upper']=resultsMAX_10
    resultOBJ[sampleOrder]['10_lower']=resultsMIN_10
######################## 35 ##############################
    file_path = './resultsMAX_35.txt'
    # Open and read the file
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # Process the lines to extract the numerical values after the '=' sign
    resultsMAX_35 = []
    for line in lines:
        if '=' in line:
            # Split the line by '=' and strip leading/trailing spaces
            value = line.split('=')[1].strip()
            # Convert the value to a float and append to the resultsMAX_35 list
            resultsMAX_35.append(float(value))

    file_path = './resultsMIN_35.txt'
    # Open and read the file
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # Process the lines to extract the numerical values after the '=' sign
    resultsMIN_35 = []
    for line in lines:
        if '=' in line:
            # Split the line by '=' and strip leading/trailing spaces
            value = line.split('=')[1].strip()
            # Convert the value to a float and append to the resultsMIN_35 list
            resultsMIN_35.append(float(value))

    resultOBJ[sampleOrder]['35_upper']=resultsMAX_35
    resultOBJ[sampleOrder]['35_lower']=resultsMIN_35
######################## 65 ##############################
    file_path = './resultsMAX_65.txt'
    # Open and read the file
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # Process the lines to extract the numerical values after the '=' sign
    resultsMAX_65 = []
    for line in lines:
        if '=' in line:
            # Split the line by '=' and strip leading/trailing spaces
            value = line.split('=')[1].strip()
            # Convert the value to a float and append to the resultsMAX_65 list
            resultsMAX_65.append(float(value))

    file_path = './resultsMIN_65.txt'
    # Open and read the file
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # Process the lines to extract the numerical values after the '=' sign
    resultsMIN_65 = []
    for line in lines:
        if '=' in line:
            # Split the line by '=' and strip leading/trailing spaces
            value = line.split('=')[1].strip()
            # Convert the value to a float and append to the resultsMIN_65 list
            resultsMIN_65.append(float(value))

    resultOBJ[sampleOrder]['65_upper']=resultsMAX_65
    resultOBJ[sampleOrder]['65_lower']=resultsMIN_65


print(resultOBJ)

with open('resultOBJ50.pkl', 'wb') as f:
    pickle.dump(resultOBJ, f)

# with open('resultOBJ500.pkl', 'wb') as f:
#     pickle.dump(resultOBJ, f)

