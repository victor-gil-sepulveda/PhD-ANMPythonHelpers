'''
Created on 25/11/2015

@author: victor
'''
import os
import pandas as pd
from nma_algo_char.common import parameter_value_to_string
import matplotlib.pyplot as plt
import seaborn as sns

if __name__ == '__main__':
    
#     folder_file = "" 
#     folders = [ line.strip() for line in open(folder_file)]
#     folders = ["Results/cc_open_steer_dispf",
#                "Results/cc_closed_rmsg_dispf",
#                "Results/cc_closed_steer_dispf",
#                "Results/cc_open_rmsg_dispf"]
# 
#     average_by = "rmsd:inc_U".split(":")
#     plot_using_all = "displacementFactor:steeringForce:rmsg".split(":")
#     
    
    folders = [
               "Results/ic_closed_iters_dispf",
               "Results/ic_open_iters_dispf"
               ]

    average_by = "rmsd:inc_U".split(":")
    plot_using_all = "displacementFactor:iters".split(":")
    
    colors = sns.color_palette("pastel", 8)
    for i,folder in enumerate(folders):
        data_path = os.path.join(folder,"data.csv")
        df = pd.read_csv(data_path, index_col=[0])
        ## locate the keys that this guy is using 
        plot_using = []
        for key in plot_using_all:
            if key in df.columns:
                plot_using.append(key)
        
        p1_keys = sorted(list(set(df[plot_using[0]])))
        p2_keys = sorted(list(set(df[plot_using[1]])))
        
        xs = []
        ys = []
        xs_std = []
        ys_std = []
                       
        labels = []        
        for p1_key in p1_keys:
            #Choose all elements with p1 key
            p1_indices = df[plot_using[0]] == p1_key
            for p2_key in p2_keys:
                p2_indices = df[plot_using[1]] == p2_key
                subdf = df.loc[p1_indices].loc[p2_indices]
                x,y = subdf.T.loc[average_by[0]].mean(), subdf.T.loc[average_by[1]].mean()
                x_err,y_err = subdf.T.loc[average_by[0]].std(), subdf.T.loc[average_by[1]].std()
                xs.append(x)
                ys.append(y)
                xs_std.append(x_err)
                ys_std.append(y_err)
                
                labels.append("%s %s"%(parameter_value_to_string(p1_key),
                                       parameter_value_to_string(p2_key)))
        plt.scatter(xs, ys, marker = 'o', c = colors[i], label = os.path.basename( folder))
        plt.errorbar(xs, ys, xerr=xs_std, yerr=ys_std, linestyle="None", c = colors[i])
        for label, x, y in zip(labels, xs, ys):
            plt.annotate(
                label, 
                xy = (x, y), xytext = (5, 5),
                textcoords = 'offset points', ha = 'right', va = 'bottom', size=6)
    lgd = plt.legend(loc='upper center', ncol=4)
    plt.show()

    # savefig bbox_extra_artists=(lgd,)
        