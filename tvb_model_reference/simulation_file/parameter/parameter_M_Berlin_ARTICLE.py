import numpy as np
import os

# states of the article : (E_L_e, E_L_i) in mV
ETATS = {'Low': (-68.0, -58.0), 'High': (-58.0, -72.0)}

class Parameter :
    def __init__(self, etat='Low'):
        path = os.path.dirname(os.path.abspath(__file__))
        E_L_e, E_L_i = ETATS[etat]
        self.etat = etat
        self.parameter_simulation={
            'path_result':'./result/synch/',
            'seed':6, 
            'save_time': 1000.0,
        }

        self.parameter_model ={
            'matteo':False,
            'delay':True, 
            'order':1,
            #parameter of the model
            'g_L':10.0,
            'E_L_e':E_L_e, #-63.0
            'E_L_i':E_L_i, #-65.0
            'C_m':200.0,
            'b_e':60.0, #4.0
            'a_e':4.0, #5.0
            'b_i':0.0,
            'a_i':0.0,
            'tau_w_e':500.0,
            'tau_w_i':1.0,
            'E_e':0.0,
            'E_i':-80.0,
            'Q_e':1.25, #1.5
            'Q_i':5.0,
            'tau_e':5.0,
            'tau_i':5.0,
            'N_tot':10000,
            'p_connect_e':0.05,
            'p_connect_i':0.05,
            'g':0.2,
            'T':5.0, #6.0
            'local_delay':3.0,
            'P_e':[-0.04983106, 0.00506355, -0.02347012, 0.00229515, -0.00041053, 0.01054705, -0.03659253, 0.00743749, 0.00126506, -0.04072161],
            'P_i':[-0.05149122, 0.00400369, -0.00835201, 0.00024142, -0.00050706, 0.00143454, -0.01468669, 0.00450271, 0.00284722, -0.0153578],
            'external_input_ex_ex':4e-3, #0.315*1e-3
            'external_input_ex_in':0.000,
            'external_input_in_ex':4e-3, #0.315*1e-3
            'external_input_in_in':0.000,
            'tau_OU':5.0,
            'weight_noise': 0.5e-3/np.sqrt(5.0), #1e-4
            'K_ext_e':400,
            'K_ext_i':0,
            #Initial condition :
            'initial_condition':{
                "E": [0.000, 0.000],"I": [0.000, 0.000],"W_e": [100.0, 100.0],"W_i": [0.0,0.0],"noise":[0.0,0.0]}
        }

        self.parameter_connection_between_region={
            ## CONNECTIVITY
            'default':False,
            'from_file':False,
            'from_folder':True,
            'from_h5':False,
            'path':path+'/../../data/QL_20120814', 
            'number_of_regions':0, #
            'tract_lengths':[],
            'weights':[],
            'speed':0.1, 
            'normalised':False
        }

        self.parameter_coupling={
            ##COUPLING
            'type':'Linear', # choice : Linear, Scaling, HyperbolicTangent, Sigmoidal, SigmoidalJansenRit, PreSigmoidal, Difference, Kuramoto
            'parameter':{'a':1.0,
                         'b':0.0}
        }

        self.parameter_integrator={
            ## INTEGRATOR
            'type':'Heun', # choice : Heun, Euler, RK4
            'stochastic': True,
            'noise_type': 'Additive', # choice : Additive
            'noise_parameter':{
                'nsig':[0.0, 0.0, 0.0, 0.0, 1.0], # first order
                'ntau':0.0,
                'dt': 0.1
                                },
            'dt': 0.1 # in ms
        }

        self.parameter_monitor= {
            'Raw':True,
            'TemporalAverage':False,
            'parameter_TemporalAverage':{
                'variables_of_interest':[0,1,2,3,4],
                'period':self.parameter_integrator['dt']*10.0
            },
            'Bold':False,
            'parameter_Bold':{
                'variables_of_interest':[0],
                'period':self.parameter_integrator['dt']*2000.0
            },
            'Ca':False,
            'parameter_Ca':{
                'variables_of_interest':[0,1,2],
                'tau_rise':0.01,
                'tau_decay':0.1
            }
        }


        self.parameter_stimulus = {
            'onset': 99.0,
            "tau": 9.0,
            "T": 99.0,
            "weights": None,
            "variables":[0]
        }
