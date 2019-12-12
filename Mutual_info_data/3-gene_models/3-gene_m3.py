# Laurent Jutras-Dubé

import numpy as np
import argparse



parser = argparse.ArgumentParser(description='Computes the mutual information between the initial and final states of the explicit model')
parser.add_argument('-p', '--param', nargs=2, type=int, help='Indices of the diffusion coefficient and noise strength')
args = parser.parse_args()



def compute_g(t):

    return min(1., np.exp(0.06*(t_osc-t)))


def dyn(var, param):

	g, s, n = param["g"], param["s"], param["n"]

	ddyn = np.zeros(var.shape)
	ddyn[0] = 1./(1.+(var[1]/s)**n)
	ddyn[1] = 1./(1.+(var[2]/s)**n)
	ddyn[2] = 1./(1.+(var[0]/s)**n)

	return (g/0.6)**5/(1.+(g/0.6)**5)*ddyn


def stat(var, param):

	g, s, n = param["g"], param["s"], param["n"]

	dstat = np.zeros(var.shape)
	dstat[0] = 1./(1.+(var[1]/s)**n) /(1.+(var[2]/s)**n)
	dstat[1] = 1./(1.+(var[2]/s)**n) /(1.+(var[0]/s)**n)
	dstat[2] = 1./(1.+(var[0]/s)**n) /(1.+(var[1]/s)**n)

	return 1./(1.+(g/0.4)**5)*dstat


def dif(var):

    dif = np.zeros(var.shape)
    d = dif_coef[dif_index]

    dif[:,:,:,0] = d*2.*(var[:,:,:,1]-var[:,:,:,0])     # Boundary
    dif[:,:,:,-1] = d*2.*(var[:,:,:,-2]-var[:,:,:,-1])  # conditions
    dif[:,:,:,1:-1] = d*(var[:,:,:,:-2]-2.*var[:,:,:,1:-1]+var[:,:,:,2:])

    return dif


def df(var, param):

    return dyn(var, param) +stat(var, param) -var +dif(var)


def compute_noise(var, param):

    noise = np.random.normal(loc=0., scale=1., size=var.shape)
    std = np.sqrt(dyn(var, param) +stat(var, param) +var)

    return std*noise



parameters = {
	"g" : 1.0,
	"s" : 0.4,
	"n" : 5.0
}


total_time = 120.
dt = 0.01
times = np.arange(0., total_time, dt)
t_osc = 10.

typ_conc = np.array([50., 100., 200., 300., 400., 500., 625., 750., 875., 1000., 1125., 1250., 1500., 1750.,
                     2000., 2500., 3000., 3500., 4000., 4500., 5000., 6000., 7500., 10000.])
noise_index = args.param[1]-1
tag_noise = True
if (noise_index < 0):    tag_noise = False

dif_coef = np.array([0., 0.01, 0.05, 0.1, 0.5])
dif_index = args.param[0]


## Sample the initial conditions

init_conc = np.genfromtxt('3-gene_m3_init_conc.txt', delimiter=',')
assert (len(init_conc[0]) == 3), "The size of the initial concentrations provided in the file init_conc.txt is not 3."

n_init = len(init_conc)
p_i = 1./n_init



## Deal with the final states

n_cell = 200
p_cell = 1./n_cell
n_run = 50
p_run = 1./n_run
rel_conc = np.zeros((n_init, 3))



## Compute the mutual information

conc = np.zeros((n_cell, n_run, n_init, 3))
conc[:] = init_conc
conc = conc.T

for t in times[1:]:

	parameters["g"] = compute_g(t)
	conc = (conc +df(conc, parameters)*dt
			+compute_noise(conc, parameters)*np.sqrt(dt/typ_conc[noise_index])*tag_noise)
	conc[conc < 0.] = 0.

conc[conc < 0.05] = 0.
total_conc = np.sum(conc, axis=0)
rel_conc += np.sum(conc/total_conc, axis=(3,2)).T


p_f_knowing_i = rel_conc*p_run*p_cell
p_f = np.sum(p_f_knowing_i*p_i, axis=0)

log_MI = np.zeros((n_init,3))
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

file = open("MI_data_m3_d"+str(dif_index)+"_n"+str(noise_index+1)+".txt", "w")
if (noise_index > -0.5):    file.write(str(typ_conc[noise_index])+","+str(mutual_info))
else:    file.write("inf,"+str(mutual_info))
file.close()
