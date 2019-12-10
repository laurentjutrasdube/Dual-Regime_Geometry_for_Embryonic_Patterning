# Laurent Jutras-Dubé

import numpy as np
import argparse



def compute_g(t):

	if (t < 4.):    g = np.ones(n_init)
	else:    g = np.zeros(n_init)

	return g


def dyn(var, param):

	g, s, n = param["g"], param["s"], param["n"]

	ddyn = np.zeros(var.shape)
	ddyn[0] = 1./(1.+(var[1]/s)**n)
	ddyn[1] = 1./(1.+(var[2]/s)**n)
	ddyn[2] = 1./(1.+(var[0]/s)**n)

	return g*g*ddyn


def stat(var, param):

	g, s, n = param["g"], param["s"], param["n"]

	dstat = np.zeros(var.shape)
	dstat[0] = 1./(1.+(var[1]/0.32)**n) /(1.+(var[2]/0.32)**n)
	dstat[1] = 1./(1.+(var[2]/0.36)**n) /(1.+(var[0]/0.36)**n)
	dstat[2] = 1./(1.+(var[0]/s)**n) /(1.+(var[1]/s)**n)

	return (1-g)*(1-g)*dstat


def df(var, param):

    return dyn(var, param) +stat(var, param) -var



parameters = {
	"g" : 1.0,
	"s" : 0.4,
	"n" : 5.0
}


total_time = 120.
dt = 0.01
times = np.arange(0., total_time, dt)



## Sample the initial conditions

init_conc = np.genfromtxt('init_conc.txt', delimiter=',')
assert (len(init_conc[0]) == 3), "The size of the initial concentrations provided in the file init_conc.txt is not 3."

n_init = len(init_conc)
p_i = 1./n_init



## Compute the mutual information

conc = init_conc.T

for t in times[1:]:

	parameters["g"] = compute_g(t)
	conc = conc +df(conc, parameters)*dt
	conc[conc < 0.] = 0.

conc[conc < 0.05] = 0.
total_conc = np.sum(conc, axis=0)
rel_conc = conc/total_conc


p_f_knowing_i = rel_conc.T
p_f = np.sum(p_f_knowing_i*p_i, axis=0)

log_MI = np.zeros((n_init, 3))
for i in range(n_init):
    for f in range(3):
        if (p_f_knowing_i[i,f] > 0.):    log_MI[i,f] = np.log2(p_f_knowing_i[i,f]/p_f[f])

mutual_info = np.sum(p_f_knowing_i*p_i*log_MI)


for k in range(len(p_f_knowing_i)):
	assert (np.sum(p_f_knowing_i, axis=1)[k] < 1.001), "Conditional probabilities do not sum to 1."
	assert (np.sum(p_f_knowing_i, axis=1)[k] > 0.999), "Conditional probabilities do not sum to 1."

assert (np.sum(p_f) < 1.001), "Probabilities of the final states do not sum to 1."
assert (np.sum(p_f) > 0.999), "Probabilities of the final states do not sum to 1."



# Write the result in an output file

file = open("MI_data_model1_stepfct.txt", "w")
file.write("inf,"+str(mutual_info))
file.close()