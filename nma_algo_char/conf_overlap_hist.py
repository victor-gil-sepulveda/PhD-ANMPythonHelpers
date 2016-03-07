import numpy
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import json
sns.set_style("whitegrid")

results = json.load(open(sys.argv[1]))

bar_width = 0.3
n_groups = len(results["cluster_ids"])
index = numpy.arange(n_groups) + bar_width/2
legend_labels = []
colors = sns.color_palette("hls", len(results["cluster_ids"]))
for j, traj in enumerate(results["id_to_path"].keys()):
    rects = plt.bar(index, 
                    results["populations"][traj], 
                    bar_width,
                    color = colors[j]
                    ) 
    index = index + bar_width
    legend_labels.append(results["id_to_path"][traj])
lgd = plt.legend(legend_labels, 
                 loc='center right', 
                 bbox_to_anchor=(1, 1))

plt.legend()
plt.show()

