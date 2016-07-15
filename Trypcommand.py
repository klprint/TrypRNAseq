# Trypcommand
import glob
import os
import sys

def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

files = glob.glob('*.txt')
scriptPath = getScriptPath()

print(files)
print(scriptPath)