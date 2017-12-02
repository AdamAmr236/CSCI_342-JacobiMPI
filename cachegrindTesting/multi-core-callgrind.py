from io import TextIOWrapper, BytesIO, StringIO
import sys
import subprocess
import re
import locale
import os

locale.setlocale( locale.LC_ALL, 'en_US.UTF-8' ) 
# parameters for jacobi iteration
threads = 4
eps = .01

# (?:elapsed time:\s*)(.+)\n
# (?:D1  misses:\s*)(\d+,?\d*)(?:\s)

# name of program and arguments, as strings

# assignment directory where callgrind files are 
path = os.getcwd()

def removeFiles():
    files = os.listdir(path)
    for file in files:
        if "callgrind.out" in file:
            os.remove(os.path.join(path,file))

def runTests(start, end, boundry=16):
    removeFiles()
    for size in range(start, end):
        cmd = ["valgrind", "--tool=callgrind", "--cache-sim=yes", "--separate-threads=yes", "../jacobi", str(size), str(size), str(threads), str(eps)]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        data = p.communicate()[0]
        string = data.decode('ascii')

        # open files
        results = {}

        files = os.listdir(path)
        for file in files:
            if "callgrind" in file:
                f = open(os.path.join(path,file), 'r')
                contents = f.read()
                f.close()
                data = re.search(r'(?:summary:\s*)(?:\d+)\s(?:\d+)\s(?:\d+)\s(?:\d+)\s(\d+)\s(\d+)\s', contents)
                if data:
                    d1misses = int(data.group(2)) + int(data.group(1))
                    results[f.name] = d1misses


        totalMisses = 0
        out = str(size) + "\t"
        for key in sorted(results.keys()):
            x = results[key]/((size*size)/threads)
            x = round(x, 3)
            out += str(x) + "\t"
            totalMisses += results[key]

        out += (str(round(totalMisses/(size*size), 3))) + "\t"

        if size * size % 16 == 0:
            out += ("yes")
        
        print(out)

        removeFiles()




runTests(10, 100)