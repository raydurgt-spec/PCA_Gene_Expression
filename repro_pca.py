import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import numpy as np


df = pd.read_csv('data/filtered.tsv.gz', sep='\t', compression='gzip', index_col=0)


labels_raw = pd.read_csv('data/class.tsv', header=None)
class_values = labels_raw[0].values



if df.shape[1] == 105:
    df = df.T

gata3_id, xbp1_id = "2625", "4404"

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

try:
    x = df[gata3_id].values
    y = df[xbp1_id].values
    
    ax1.scatter(x[class_values==0], y[class_values==0], c='black', marker='s', s=25, label='ER-')
    ax1.scatter(x[class_values==1], y[class_values==1], c='red', marker='s', s=25, label='ER+')
    ax1.set_xlabel(f"GATA3 ({gata3_id})")
    ax1.set_ylabel(f"XBP1 ({xbp1_id})")
except KeyError:
 
    print("Specific IDs not found, using first two available columns for 1a.")
    ax1.scatter(df.iloc[:, 0], df.iloc[:, 1], c=['red' if c==1 else 'black' for c in class_values], marker='s')

ax1.set_title("Figure 1a: Gene Expression")


pca = PCA(n_components=1)
pc1 = pca.fit_transform(df).flatten()


ax2.scatter(pc1, np.full_like(pc1, 2), c=['red' if c==1 else 'black' for c in class_values], s=15)
ax2.scatter(pc1[class_values==0], np.full_like(pc1[class_values==0], 1), c='black', s=15)
ax2.scatter(pc1[class_values==1], np.full_like(pc1[class_values==1], 0), c='red', s=15)

ax2.set_yticks([2, 1, 0])
ax2.set_yticklabels(['All', 'ER-', 'ER+'])
ax2.set_xlabel("Projection onto PC1")
ax2.set_title("Figure 1c: PC1 Projection")

plt.tight_layout()
plt.savefig('repro_figure_1.png')
print("--- Success! Created repro_figure_1.png ---")