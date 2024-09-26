# import os

# # 设置 GAMS 安装目录
# os.environ['GAMSDIR'] = '/mnt/d/gams/47'

# import gdxpds
# gdxpds.load_gdxcc(gams_dir=os.environ['GAMSDIR'])

# 现在导入其他库
import pandas as pd
import subprocess
import concurrent.futures

def run_gams(gams_file):
    print(f"开始运行: {gams_file}")  # 打印任务开始信息
    result = subprocess.run(['/mnt/d/gams/47/gams.exe', gams_file], capture_output=True, text=True)
    print(f"完成运行: {gams_file}")  # 打印任务结束信息
    return result.stdout


gams_files = ['LULinearPortfolioShortLoopMax.gms', 'LULinearPortfolioShortLoopMin.gms']

with concurrent.futures.ProcessPoolExecutor() as executor:
    # 提交多个任务并行运行
    results = list(executor.map(run_gams, gams_files))

file_path = './resultsMAX_1.txt'

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