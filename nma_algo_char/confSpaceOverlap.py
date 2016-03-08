"""
Created on Mar 7, 2016

@author: victor
"""
import sys
import json
import math
import numpy
from scipy.stats import entropy
from numpy.linalg import norm
import anmichelpers.tools.tools as tools
from pyproct.data.handler.protein.proteinEnsembleDataLoader import ProteinEnsembleDataLoader
from pyproct.data.matrix.protein.cases.rmsd.cartesiansCase import RMSDMatrixBuilder
from pyproct.clustering.algorithms.kmedoids.kMedoidsAlgorithm import KMedoidsAlgorithm
from pyproct.data.handler.dataSource import DataSource
from pyproct.driver.parameters import ProtocolParameters

def get_matrix_params(prot_type, fit_selection):
    if prot_type == "src":
        return ProtocolParameters({
            "matrix": {
                "method": "rmsd::ensemble",
                "parameters": {
                    "calculator_type": "QCP_OMP_CALCULATOR",
                    "calculator_options":{
                        "number_of_threads": 6
                    },
                    "fit_selection": "(name CA) and ((resnum 266 to 324) or (resnum 336 to 345) or (resnum 363 to 406) or (resnum 433 to 454) or (resnum 481 to 524))"
                }
            }})
        
    if prot_type == "ubi":
        return ProtocolParameters({
        "matrix": {
            "method": "rmsd::ensemble",
            "parameters": {
                "calculator_type": "QCP_OMP_CALCULATOR",
                "calculator_options":{
                    "number_of_threads":6
                },
                "fit_selection": "(name CA) and (resnum 1 to 70)"
            }
        }
    })
    else:
        return None

def calc_interdomain_distances(pdb, elements): 
    res_coords = pdb.select("(name CA) and (resid 277 or resid 387)").getCoordsets()[elements]
    ds = []
    for cys_coords, leu_coords in res_coords:
        ds.append(tools.norm(leu_coords-cys_coords))
    return numpy.mean(ds), numpy.std(ds)

def get_matrix(prot_type, sources, fit_selection):
    loader = ProteinEnsembleDataLoader(get_matrix_params(prot_type, fit_selection))
    
    e_ranges = []
    for source in sources:
        e_range = loader.load(source)
        e_ranges.append(e_range)
    ensemble_data = loader.close()
    
    class mock_handler:
        def __init__(self, data):
            self.inner_data = data
        def get_data(self):
            return self.inner_data
    
    matrix = RMSDMatrixBuilder.build(mock_handler(ensemble_data), get_matrix_params(prot_type, fit_selection))
    
    return ensemble_data, matrix, e_ranges

def process_ensemble( prot_type, representatives, clusters_per_representative,  element_range, the_matrix):
    min_distances_per_rep = {}
    ens_populations_per_rep = {}
    for r in representatives:
        ens_populations_per_rep[r] = 0
    
    ens_cluster_per_rep = {}
    for element in element_range:
        # Get the closer representative
        distances = []
        for r in representatives:
            distances.append((the_matrix[r,element], r))
        distances.sort()
        _, closer_r = distances[0]
        
        ens_populations_per_rep[closer_r] += 1
        
        try:
            ens_cluster_per_rep[closer_r].append(element)
        except:
            ens_cluster_per_rep[closer_r] = [element]
        
        # Get the minimum dist of these element with the elements of ref. cluster
        distances = []
        c = clusters_per_representative[closer_r]
        for c_el in c.all_elements:
            distances.append(the_matrix[element, c_el])
        try:
            min_distances_per_rep[closer_r].append(min(distances))
        except:
            min_distances_per_rep[closer_r] = [min(distances)]
    
    return ens_populations_per_rep, min_distances_per_rep, ens_cluster_per_rep

def smoothed(distribution, small_value = 1.0e-8):
    """
    Applies a smoothing process to the distribution.
    See http://mathoverflow.net/questions/72668/how-to-compute-kl-divergence-when-pmf-contains-0s
    for an explanation about the problem and the solution.
     
    @param distribution: distribution to be smoothed
    @param small_value: value to be set to those bins with 0 probability
     
    @return: The smoothed distribution.
    """
    total_number_of_samples = len(distribution)
    samples_in_distrib = numpy.count_nonzero(distribution)
    if samples_in_distrib > 0:
        pc = small_value * (total_number_of_samples - samples_in_distrib) / samples_in_distrib
        smoothed_distrib = numpy.empty(len(distribution))
        for i in range(len(distribution)):
            if distribution[i] == 0:
                smoothed_distrib[i] = small_value
            else:
                smoothed_distrib[i] = distribution[i] - pc
        return numpy.array(smoothed_distrib)
    else:
        return distribution

def JSD(P, Q):
    """
    Calculates the Jensen-Shannon divergence as a metric (sq_root)
    See: http://www.researchgate.net/publication/3084774_A_new_metric_for_probability_distributions
    """
    _P = P / norm(P, ord=1)
    _Q = Q / norm(Q, ord=1)
    _M = 0.5 * (_P + _Q)
    return math.sqrt(0.5 * (entropy(_P, _M) + entropy(_Q, _M)))
    
def get_avg_distances(min_distances):
    avgs = {}
    all_distances = []
    for rep in min_distances:
        avgs[rep] =  (numpy.mean(min_distances[rep]),
                     numpy.std(min_distances[rep]))
        all_distances.extend(min_distances[rep])
    return avgs, (numpy.mean(all_distances),numpy.std(all_distances))

if __name__ == '__main__':

    prot_type = sys.argv[1]

    if prot_type == "src":
#         sources = {
#                    "MD": {
#                          "source": "src_md_ca.pdb",
#                          "base_selection": "not (name OXT) and heavy"
#                    },
#                    "IC": {
#                          "source": "src_ic_ca.pdb",
#                          "base_selection":"not (name OXT) and heavy"
#                    },
#                    "CC": {
#                          "source": "src_cc_ca.pdb",
#                          "base_selection":"not (name OXT) and heavy"
#                    }
#         }
        sources = {
                   "MD": {
                         "source": "src_md_compressed.pdb",
                         "base_selection": "not (name OXT) and heavy"
                   },
                   "IC": {
                         "source": "src_ic_compressed.pdb",
                         "base_selection":"not (name OXT) and heavy"
                   },
                   "CC": {
                         "source": "cc_src_compressed.pdb",
                         "base_selection":"not (name OXT) and heavy"
                   }
        }
        fit_selection = "(name CA) and ((resnum 266 to 324) or (resnum 336 to 345) or (resnum 363 to 406) or (resnum 433 to 454) or (resnum 481 to 524))"
    
    if prot_type == "ubi":
        sources = {
                   "MD":{
                        "source":"ubi_md_compressed.pdb",
                        "base_selection": "not (name OXT) and heavy and (resnum 1 to 73)"
                    },
                   "IC":{
                            "source":"ic_ubi_compressed.pdb" ,
                            "base_selection": "not (name OXT) and heavy (resnum 1 to 73)"
                    },
                   "CC":{
                            "source": "cc_ubi_compressed.pdb",
                            "base_selection": "not (name OXT) and heavy (resnum 1 to 73)"
                    }
          
        }
        fit_selection = "(name CA) and (resnum 1 to 70)"
    
    
    results = {}
    ref_ensemble, ref_matrix, _ = get_matrix(prot_type, [DataSource(sources["MD"])], fit_selection)
    first_merged_ensemble, matrix1, e_ranges_ref_ens1 = get_matrix(prot_type, [DataSource(sources["MD"]), DataSource(sources["IC"])], fit_selection)
    second_merged_ensemble, matrix2, e_ranges_ref_ens2 = get_matrix(prot_type, [DataSource(sources["MD"]), DataSource(sources["CC"])], fit_selection)
    for k in range(1,20) + range(20,100,5) + range(100,200, 10) + range(200, 500, 50)+range(500, 1001, 100):
        clustering = KMedoidsAlgorithm(ref_matrix).perform_clustering({"k": k, "seeding_type":"EQUIDISTANT"})
        print "Clustering performed for k = ",k
        
        representatives = [c.prototype for c in  clustering.clusters]
        ref_populations = dict(zip(representatives,[len(c.all_elements) for c in  clustering.clusters]))
        clusters_per_representative = dict(zip(representatives,clustering.clusters))
        if prot_type == "src":
            ref_interdomain_dists = []
            for r in representatives:
                dists = calc_interdomain_distances(ref_ensemble.get_all_elements(), clusters_per_representative[r].all_elements)
                ref_interdomain_dists.append(dists)
        
        print "\tProcessing IC"
        ens1_populations, ens1_min_distances, ens1_cluster_per_rep = process_ensemble(prot_type, representatives, clusters_per_representative, e_ranges_ref_ens1[1], matrix1)
        if prot_type == "src":
            ens1_interdomain_dists = []
            for r in representatives:
                if r in ens1_cluster_per_rep.keys():
                    dists = calc_interdomain_distances(first_merged_ensemble.get_all_elements(), ens1_cluster_per_rep[r])
                    ens1_interdomain_dists.append(dists)
                else:
                    ens1_interdomain_dists.append((0.,0.))
                
        print "\tProcessing CC"    
        ens2_populations, ens2_min_distances, ens2_cluster_per_rep = process_ensemble(prot_type, representatives, clusters_per_representative, e_ranges_ref_ens2[1], matrix2)   
        if prot_type == "src":
            ens2_interdomain_dists = []
            for r in representatives:
                if r in ens2_cluster_per_rep.keys():
                    dists = calc_interdomain_distances(second_merged_ensemble.get_all_elements(), ens2_cluster_per_rep[r])
                    ens2_interdomain_dists.append(dists)
                else:
                    ens2_interdomain_dists.append((0.,0.))
        
        results[k] = {
                      "MD_pop": list(numpy.array([ref_populations[r]  for r in representatives])/float(len(e_ranges_ref_ens1[0]))),
                      "IC_pop": list(numpy.array([ens1_populations[r] for r in representatives])/float(len(e_ranges_ref_ens1[1]))),
                      "CC_pop": list(numpy.array([ens2_populations[r] for r in representatives])/float(len(e_ranges_ref_ens2[1]))),
                      "representatives":representatives
                      }
        if prot_type == "src":
            results[k]["MD_dists"] = ref_interdomain_dists
            results[k]["IC_dists"] = ens1_interdomain_dists
            results[k]["CC_dists"] = ens2_interdomain_dists
        
        ens1_avgs, ens1_global_avg = get_avg_distances(ens1_min_distances)
        results[k]["IC_dist_avgs"] = ens1_avgs
        results[k]["IC_dist_global_avg"] = ens1_global_avg
        
        ens2_avgs, ens2_global_avg = get_avg_distances(ens2_min_distances)
        results[k]["CC_dist_avgs"] = ens2_avgs
        results[k]["CC_dist_global_avg"] = ens2_global_avg
        
        results[k]["JSD"] = {
                             "md_ic": JSD(smoothed(results[k]["MD_pop"]), smoothed(results[k]["IC_pop"])),
                             "md_cc": JSD(smoothed(results[k]["MD_pop"]), smoothed(results[k]["CC_pop"])),
                             "ic_cc": JSD(smoothed(results[k]["IC_pop"]), smoothed(results[k]["CC_pop"]))
                             }
        
    open("SRC_conf_overlap_result.json", "w").write(
            json.dumps(results, sort_keys = False, indent = 4, separators = (',', ': '))
        )
    
    