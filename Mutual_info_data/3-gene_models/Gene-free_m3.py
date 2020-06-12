# Laurent Jutras-Dubé

import numpy as np
import argparse



parser = argparse.ArgumentParser(description='Computes the mutual information between the initial and final states of the toy model')
parser.add_argument('-p', '--param', nargs=1, type=int, help='Index of the noise strength')
args = parser.parse_args()



## Functions for the differential equations defining the D/S system

def compute_g(t):

	par = np.exp(0.036*(t_osc-t))
	if (par > 1.):    par = 1.
	elif (par < 0.):    par = 0.

	return par


def stat(var):

	dstat = np.zeros(vec.shape)
	dstat[0] = var[0]*(1.-var[0]*var[0])
	dstat[1] = -var[1]

	return dstat


def dyn(var):

	r = np.sqrt(var[0]**2 +var[1]**2)

	ddyn = np.zeros(vec.shape)
	ddyn[0] = var[0]*r*(1.-r) -var[1]
	ddyn[1] = var[1]*r*(1.-r) +var[0]

	return ddyn


def df(var, par):

	return par**3*dyn(var) +(1.-par)**3*stat(var) -par*(1.-par)*var




init_time = 0.
total_time = 150.
dt = 0.01
times = np.arange(init_time, init_time+total_time, dt)
t_osc = 20.

typ_conc = np.array([50., 100., 200., 300., 400., 500., 625., 750., 875., 1000., 1125., 1250., 1500., 1750.,
                     2000., 2500., 3000., 3500., 4000., 4500., 5000., 6000., 7500., 10000.])
noise_index = args.param[0]-1
tag_noise = True
if (noise_index < 0):    tag_noise = False




## Sample the initial conditions

n_init = 50
p_i = 1/n_init
init_vec = np.zeros((n_init, 2))

dphi = 2*np.pi*p_i
for i in range(n_init):
	init_vec[i,0] = np.cos(i*dphi)
	init_vec[i,1] = np.sin(i*dphi)




## Deal with the final states

n_cell = 1000
p_cell = 1/n_cell
x_f = np.zeros((2, n_init))




## Compute the mutual information

vec = np.zeros((n_cell, n_init, 2))
vec[:] = init_vec
vec = vec.T



# results = [vec]



for t in times[1:]:

	par = compute_g(t)
	noise = np.random.normal(loc=0., scale=1., size=vec.shape)
	vec = vec +df(vec, par)*dt +noise*np.sqrt(dt/typ_conc[noise_index])*tag_noise



	# results.append(vec)



vec[0][vec[0] < 0.] = 0.
vec[0][vec[0] > 0.] = 1.
x_f[0] = np.sum(vec[0], axis=1)
x_f[1] = n_cell -x_f[0]
x_f = x_f.T


p_f_knowing_i = x_f*p_cell
p_f = np.sum(p_f_knowing_i*p_i, axis=0)

log_MI = np.zeros((n_init, 2))
for i in range(n_init):
	for f in range(2):
		if (p_f_knowing_i[i,f] > 0):    log_MI[i,f] = np.log2(p_f_knowing_i[i,f]/p_f[f])


mutual_info = np.sum(p_f_knowing_i*p_i*log_MI)


for k in range(len(p_f_knowing_i)):
	assert (np.sum(p_f_knowing_i, axis=1)[k] < 1.001), "Conditional probabilities do not sum to 1."
	assert (np.sum(p_f_knowing_i, axis=1)[k] > 0.999), "Conditional probabilities do not sum to 1."

assert (np.sum(p_f) < 1.001), "Probabilities of the final states do not sum to 1."
assert (np.sum(p_f) > 0.999), "Probabilities of the final states do not sum to 1."




## Write the result in an output file

file = open("MI_data_gene-free_m3_"+str(noise_index+1)+".txt", "w")
if (noise_index > -0.5):    file.write(str(typ_conc[noise_index])+","+str(mutual_info))
else:    file.write("inf,"+str(mutual_info))
file.close()