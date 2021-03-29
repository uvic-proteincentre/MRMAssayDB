from collections import Counter
from itertools import combinations
import ast
from operator import itemgetter
import operator
import pandas as pd
import csv
import numpy as np

import matplotlib.pyplot as plt
import venn


labels = venn.get_labels([range(10), range(5, 15), range(3, 8), range(8, 17)], fill=['number', 'logic'])
fig, ax = venn.venn4(labels, names=['list 1', 'list 2', 'list 3', 'list 4'])
fig.show()