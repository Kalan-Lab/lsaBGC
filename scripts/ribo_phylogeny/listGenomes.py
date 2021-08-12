import os
import sys

for f in os.listdir('../Proteomes/'):
    if f.endswith('.faa'):
        print(f.split('.faa')[0])
