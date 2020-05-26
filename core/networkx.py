import networkx
import pandas as pd
import sys
import glob
import os
import logging
import seaborn as sns
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patches as patches
from matplotlib.legend_handler import HandlerPatch
import math
from collections import OrderedDict
import itertools
pd.options.mode.chained_assignment = None
import datetime


path="/Users/denizdogan/Desktop/xx/MINNIE_run_v2 2/all/data/05_common_intra-inter/25_freq_filtered/25_freq/inter/protein-dna"
df=pd.read_csv(f'{path}/Sox4_hbond_inter.csv')


df['donor_acceptor'].str.contains('6_C')
