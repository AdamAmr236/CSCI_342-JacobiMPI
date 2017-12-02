from io import TextIOWrapper, BytesIO, StringIO
import sys
import subprocess
import re
import locale
import os

locale.setlocale( locale.LC_ALL, 'en_US.UTF-8' )
# parameters for jacobi iteration
threads = 4
eps = .0005

gridSizes = {200, 400, 800, 1000, 1500, 2000, 3000, 4000, 5000, 15000}

cycleResults = {}
clockResults = {}

for i in range(15):
    for size in sorted(gridSizes):
        cmd = ["../jacobi", str(size), str(size), str(threads), str(eps)]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        data = p.communicate()[0]
        string = data.decode('ascii')

        cycleRegex = r'(?:number of cycles:\s*)(\d+)'
        cycleData = re.search(cycleRegex, string)

        clockRegex = r'(?:wall clock time:\s*)(.*)'
        clockData = re.search(clockRegex, string)

        if clockResults.get(size, float("inf")) > float(clockData.group(1)):
            clockResults[size] = float(clockData.group(1))

        if cycleResults.get(size, float("inf")) > int(cycleData.group(1)):
            cycleResults[size] = int(cycleData.group(1))



    print("iteration: " + str(i))


for key in sorted(cycleResults.keys()):
    print(str(key) + "\t" + str(cycleResults[key]) + "\t" +  str(clockResults[key]))
