import numpy as np
import EM_Files as newfiles
import EM_files as oldfiles
import EM_Algorithm as newEM
import EM_binary as oldEM
import PrettyPlots as pplots
import EM_Plots as Plots
import Simulate as sim
from EM_lib import real_cluster_probs_to_obs
from calculate_loglik import calculate_ll
import EM_Functions
from importlib import reload
import loglikelihood as ll
from importlib import reload


mu, pi = sim.LoadSimInfo("/Users/vcorbin/Workspace/RNA_structure/Method_paper/OC57_sim_95_5_BV_clusters_info.txt")
piobs = real_cluster_probs_to_obs(pi, mu, normalise=True)
muzero = mu.copy()
mu_sum = mu[0]*piobs[0] + mu[1]*piobs[1]
mu_sum = np.array([mu_sum, mu_sum])
muzero[mu_sum < 0.005] = 0.0

Xzero = newfiles.Load_BitVectors("/Users/vcorbin/Workspace/RNA_structure/Method_paper/OC57_sim_95_5_BV.txt", 0.05, 0.005, "YES", "/Users/vcorbin/Workspace/RNA_structure/Method_paper/simulate/", nmuts_thresh=10)
log_like_list_zero1, final_mu_zero1, final_obs_pi_zero1, final_real_pi_zero1, resps_zero1, BIC_zero1 = newEM.Run_EM(Xzero, 1, 10, 0.05)
log_like_list_zero2, final_mu_zero2, final_obs_pi_zero2, final_real_pi_zero2, resps_zero2, BIC_zero2 = newEM.Run_EM(Xzero, 2, 1000, 0.05)
ll_zero2 = ll.ll(Xzero.BV_Matrix, Xzero.BV_Abundance, final_mu_zero2, final_obs_pi_zero2)
llsim_zero = ll.ll(Xzero.BV_Matrix, Xzero.BV_Abundance, muzero, piobs)

calc_inds = EM_Functions.calc_matrixIndices(Xzero.BV_Matrix.shape[0], 2,18)
resps_sim, log_like_sim, denom_sim = newEM.Exp_Step(Xzero, 2, mu, piobs, calc_inds)

Xnew = newfiles.Load_BitVectors("/Users/vcorbin/Workspace/RNA_structure/Method_paper/OC57_sim_95_5_BV.txt", 0.05, 0.0, "YES", "/Users/vcorbin/Workspace/RNA_structure/Method_paper/simulate/", nmuts_thresh=10)
log_like_list_new1, final_mu_new1, final_obs_pi_new1, final_real_pi_new1, resps_new1, BIC_new1 = newEM.Run_EM(Xnew, 1, 10, 0.05)
log_like_list_new2, final_mu_new2, final_obs_pi_new2, final_real_pi_new2, resps_new2, BIC_new2 = newEM.Run_EM(Xnew, 2, 1000, 0.05)
ll_new2 = ll.ll(Xnew.BV_Matrix, Xnew.BV_Abundance, final_mu_new2, final_obs_pi_new2)
llsim_new = ll.ll(Xnew.BV_Matrix, Xnew.BV_Abundance, mu, piobs)


Xold, labels = oldfiles.load_binary("/Users/vcorbin/Workspace/RNA_structure/Method_paper/OC57_sim_95_5_BV.txt", max_mut_prob=10.0)
lls_old1, mu_old1, cps_old1 = oldEM.runEM(Xold, 1, labels, threading=True, max_its=10)
BIC_old1 = oldfiles.bic(lls_old1[-1], 1, len(mu_old1[0]), Xold.N.sum())
lls_old2, mu_old2, cps_old2 = oldEM.runEM(Xold, 2, labels, threading=True, max_its=1000)
BIC_old2 = oldfiles.bic(lls_old2[-1], 2, len(mu_old2[0]), Xold.N.sum())

pi_obs_old2 = real_cluster_probs_to_obs(cps_old2[-1], mu_old2, normalise=True)
ll_old2 = ll.ll(Xold.bitMatrix.values, Xold.N, mu_old2, pi_obs_old2)
llsim_old = ll.ll(Xold.bitMatrix.values, Xold.N, mu, piobs)

loglik, pi_real = calculate_ll(Xold, mu, obs_cluster_probs=piobs, return_cluster_probs=True)
resp2, loglik2 = oldEM.Expectation_step(Xold, mu, np.log(pi))



-----------------------------------------------------dmuts=0--------------------------------------------
Xd0, labels = oldfiles.load_binary("/Users/vcorbin/Workspace/RNA_structure/Method_paper/OC57_sim_95_5_BV.txt", max_mut_prob=10.0, dmuts=0)
lls_d0_1, mu_d0_1, cps_d0_1 = oldEM.runEM(Xd0, 1, labels, threading=True, max_its=10, dmuts=0)
BIC_d0_1 = oldfiles.bic(lls_d0_1[-1], 1, len(mu_d0_1[0]), Xd0.N.sum())
lls_d0_2, mu_d0_2, cps_d0_2 = oldEM.runEM(Xd0, 2, labels, threading=True, max_its=1000, dmuts=0)
BIC_d0_2 = oldfiles.bic(lls_d0_2[-1], 2, len(mu_d0_2[0]), Xd0.N.sum())

pi_obs_d0_2 = real_cluster_probs_to_obs(cps_d0_2[-1], mu_d0_2, dmuts=0, normalise=True)
ll_d0_2 = ll.ll(Xd0.bitMatrix.values, Xd0.N, mu_d0_2, pi_obs_d0_2, dmuts=0)
llsim_d0 = ll.ll(Xd0.bitMatrix.values, Xd0.N, mu, piobs, dmuts=0)

