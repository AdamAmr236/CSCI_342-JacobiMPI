from io import TextIOWrapper, BytesIO, StringIO
import sys
import subprocess
import re
import locale
import os

locale.setlocale( locale.LC_ALL, 'en_US.UTF-8' ) 
# parameters for jacobi iteration
threads = 1
eps = .01

# (?:elapsed time:\s*)(.+)\n
# (?:D1  misses:\s*)(\d+,?\d*)(?:\s)

# name of program and arguments, as strings

# assignment directory where callgrind files are 
path = os.getcwd()

def removeFiles():
    files = os.listdir(path)
    for file in files:
        if "cachegrind.out" in file:
            os.remove(os.path.join(path,file))

removeFiles()

for size in range(5, 200):
    cmd = ["valgrind", "--tool=cachegrind", "../jacobi", str(size), str(size), str(threads), str(eps)]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    data = p.communicate()[0]
    string = data.decode('ascii')

    removeFiles()


    d1missesStr = re.search(r'(?:D1  misses:\s*)(\d+,?\d*,?\d*,?\d*)(?:\s)', string).group(1)
    d1missesInt = locale.atoi(d1missesStr)
    d3missesStr = re.search(r'(?:LL misses:\s*)(\d+,?\d*,?\d*,?\d*)(?:\s)', string).group(1)
    d3missesInt = locale.atoi(d3missesStr)
    print(str(size) + "\t" + str(size * size) + "\t" + str(round(d1missesInt/(size*size), 3)))
    #print(str(size) + "\t" + d1missesStr + "\t" + str(d1missesInt/(size*size)) + "\t")