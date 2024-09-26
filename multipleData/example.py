# import os

# # 设置 GAMS 安装目录
# os.environ['GAMSDIR'] = '/mnt/d/gams/47'

# import gdxpds
# gdxpds.load_gdxcc(gams_dir=os.environ['GAMSDIR'])

# 现在导入其他库
import pandas as pd
import subprocess
#result = subprocess.run(['/mnt/d/gams/47/gams.exe', 'model.gms'])
result = subprocess.run(['/mnt/d/gams/47/gams.exe', 'LULinearPortfolioShortLoopMax.gms'])

file_path = './results.txt'

# Open and read the file
with open(file_path, 'r') as file:
    lines = file.readlines()

# Process the lines to extract the numerical values after the '=' sign
results = []
for line in lines:
    if '=' in line:
        # Split the line by '=' and strip leading/trailing spaces
        value = line.split('=')[1].strip()
        # Convert the value to a float and append to the results list
        results.append(float(value))

print(results)