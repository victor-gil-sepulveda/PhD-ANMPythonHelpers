'''
Created on Mar 8, 2016

@author: victor
'''
import json
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import numpy

sns.set_style("whitegrid")

if __name__ == '__main__':
    results = json.load(open(sys.argv[1]))
    
    # get JSD values per k
    ic_to_cc = []
    md_to_ic = []
    md_to_cc = []
    
    for str_k in results:
        k = int(str_k)
        ic_to_cc.append((k,results[str_k]["JSD"]["ic_cc"]))
        md_to_ic.append((k,results[str_k]["JSD"]["md_ic"]))
        md_to_cc.append((k,results[str_k]["JSD"]["md_cc"]))

    ic_to_cc.sort()
    md_to_ic.sort()
    md_to_cc.sort()

#     ax1 = plt.subplot2grid((4,1), (0,0), rowspan=3)
#     
#     for stuff_to_plot in [("JSD(MD, IC)", md_to_ic),("JSD(MD, CC)", md_to_cc),("JSD(IC, CC)", ic_to_cc)]:
#         label, data = stuff_to_plot[0], numpy.array(stuff_to_plot[1])
#         ax1.plot(data.T[0], data.T[1], label = label,  marker = "o")
#     ax1.set_xscale("log", nonposx='clip')
#     ax1.get_xaxis().set_visible(False)
#     plt.legend()
#     
#     # get values for the rmsd average etc
#     ic_rmsd = []
#     cc_rmsd = []
#     
#     for str_k in results:
#         k = int(str_k)
#         ic_rmsd.append((k, results[str_k]["IC_dist_global_avg"][0],results[str_k]["IC_dist_global_avg"][1]))
#         cc_rmsd.append((k,results[str_k]["CC_dist_global_avg"][0],results[str_k]["CC_dist_global_avg"][1]))
#     ic_rmsd.sort()
#     cc_rmsd.sort()
#     
#     
#     ax2 = plt.subplot2grid((4,1), (3,0), colspan=1)
# 
#     for stuff_to_plot in [("IC avg. RMSD", ic_rmsd),("CC avg. RMSD", cc_rmsd)]:
#         label, data = stuff_to_plot[0], numpy.array(stuff_to_plot[1])
#         ax2.errorbar(data.T[0], data.T[1], yerr = data.T[2], label = label,  marker = "o")
#     plt.gca().set_xscale("log", nonposx='clip')
#     plt.legend()
#     plt.show()
#     plt.close()
    
    chosen_k = int(sys.argv[2])
    
    plot_distances = "IC_dists" in results[str(chosen_k)]
    
    if not plot_distances:
        bar_width = 0.3
        n_groups = len(results[str(chosen_k)]["MD_pop"])
        index = numpy.arange(n_groups) + bar_width/2
        legend_labels = []
        colors = sns.color_palette("hls", 3)
        for j, ens_label in enumerate(["MD_pop", "IC_pop", "CC_pop"]):
            print results[str(chosen_k)][ens_label]
            rects = plt.bar(index, 
                            results[str(chosen_k)][ens_label], 
                            bar_width,
                            color = colors[j]
                            ) 
            index = index + bar_width
            legend_labels.append(ens_label)
        lgd = plt.legend(legend_labels, 
                         loc='center right', 
                         bbox_to_anchor=(1, 1))
        
        plt.legend()
        plt.show()
    else:
        ax1 = plt.subplot2grid((4,1), (0,0), rowspan=3)
        bar_width = 0.3
        n_groups = len(results[str(chosen_k)]["MD_pop"])
        index = numpy.arange(n_groups) - 0.5
        legend_labels = []
        colors = sns.color_palette("hls", 3)
        for j, ens_label in enumerate(["MD_pop", "IC_pop", "CC_pop"]):
            print results[str(chosen_k)][ens_label]
            rects = ax1.bar(index, 
                            results[str(chosen_k)][ens_label], 
                            bar_width,
                            color = colors[j]
                            ) 
            index = index + bar_width
            legend_labels.append(ens_label)
        ax1.set_xlim(xmin=-0.5, xmax=chosen_k+0.5)
        plt.legend(legend_labels, 
                         loc='center right', 
                         bbox_to_anchor=(1, 1))
        
        ax2 = plt.subplot2grid((4,1), (3,0), colspan=1)
        for label, data, color in [("MD interd. dist.", results[str(chosen_k)]["MD_dists"], colors[0]),
                                   ("IC interd. dist.", results[str(chosen_k)]["IC_dists"], colors[1]),
                                   ("CC interd. dist.", results[str(chosen_k)]["CC_dists"], colors[2])]:
            data = numpy.array(data)
            nonzeros = numpy.nonzero(data.T[0])[0]
            print nonzeros, label, data
            ax2.errorbar(x = nonzeros, 
                         y = data.T[0][nonzeros], 
                         yerr = data.T[1][nonzeros], 
                         label = label, 
                         fmt='o', 
                         color = color)
        ax2.set_xlim(xmin=-0.5, xmax=chosen_k+0.5)
        plt.legend()
        plt.show()
        plt.close()
            