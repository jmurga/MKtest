from adapter import extensions
from adapter import importData
import pandas as pd
import numpy as np


rawData = importData.MultiFasta('/home/jmurga/subset.fa','standard',None,50)

m = rawData.sequencesToMatrix()
s = rawData.uSfsFromFasta(m)
sfs, div =rawData.formatSfs(s,m)


analysis = extensions.Mk(sfs,div)
print(analysis.createTableFull())
