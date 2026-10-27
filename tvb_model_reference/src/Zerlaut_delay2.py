import matplotlib.pyplot as plt
import numpy
from tvb.simulator.lab import *
from tvb.simulator.history import SparseHistory
from tvb_model_reference.src.Zerlaut import Zerlaut_adaptation_first_order as Model1, Zerlaut_adaptation_second_order as Model2
from tvb.basic.neotraits.api import NArray, Range, Final, List


class MyHistory(SparseHistory): 
    dt = None
    step = -1
    def query_sparse(self, step):
        self.step = step 
        return super().query_sparse(step)
    def get_delayed_state(self, tau):
        tau_steps = int(tau / self.dt) 
        assert tau_steps < self.n_time  
        return self.buffer[ (self.step - 1 - tau_steps) % self.n_time ] 

class MySim(simulator.Simulator): 
    def configure(self):
        super().configure()
        self.model.history = self.history
        self.history.dt = self.integrator.dt
        return self
    def _configure_history(self, initial_conditions=None):
        self.history = MyHistory.from_simulator(self, initial_conditions) 
        

class Zerlaut_adaptation_first_order(Model1): 

    local_delay =  NArray(
        label=":math:`\ntau local`",
        default=numpy.array([5]),
        domain=Range(lo=0., hi=50.0, step=1.0),
        doc="""local tau""")

    cvar = numpy.array([0, 1], dtype=int)
    
    epsilon = 0.05 
    stim_kHz = None 
    stim_noeud = None 

    def dfun(self, state_variables, coupling, local_coupling=0.00): 

        E, I, W_e, W_i, noise = state_variables
        tau = float(numpy.squeeze(self.local_delay))
        delayed_long = self.history.get_delayed_state((1 + self.epsilon) * tau)
        delayed_short = self.history.get_delayed_state((1 - self.epsilon) * tau)
        E_long, I_long = delayed_long[0], delayed_long[1]
        E_short, I_short = delayed_short[0], delayed_short[1]
        derivative = numpy.zeros_like(state_variables)

        c_0 = coupling[0, :]

        Fe_ext = c_0 + self.weight_noise * noise
        if self.stim_kHz is not None: 
            k = self.history.step
            if 0 <= k < len(self.stim_kHz):
                Fe_ext = Fe_ext.copy()
                Fe_ext[self.stim_noeud] += self.stim_kHz[k]
        input_e = Fe_ext + self.external_input_ex_ex
        input_i = Fe_ext + self.external_input_in_ex
        input_e = numpy.where(input_e < 0, 0.0, input_e) 
        input_i = numpy.where(input_i < 0, 0.0, input_i)
        Fi_ext = 0.0 * Fe_ext

        derivative[0] = (self.TF_excitatory(
            E_long, I_short, input_e,
            Fi_ext + self.external_input_ex_in,
            W_e)-E)/self.T

        derivative[1] = (self.TF_inhibitory(
            E_short, I_long, input_i,
            Fi_ext + self.external_input_in_in,
            W_i)-I)/self.T
        
        # Adaptation excitatory
        zero = 0.0 * Fe_ext
        mu_V, sigma_V, T_V = self.get_fluct_regime_vars(
                E, I, zero, zero, W_e, self.Q_e, self.tau_e, self.E_e,
                self.Q_i, self.tau_i, self.E_i,
                self.g_L, self.C_m, self.E_L_e, self.N_tot, self.p_connect_e,
                self.p_connect_i, self.g, self.K_ext_e, self.K_ext_i)
        derivative[2] = -W_e/self.tau_w_e+self.b_e*E+self.a_e*(mu_V-self.E_L_e)/self.tau_w_e
        # Adaptation inhibitory
        derivative[3] = 0.0 * W_i
        derivative[4] = -noise/self.tau_OU

        return derivative


class Zerlaut_adaptation_second_order(Model2): 

    local_delay =  NArray(
        label=":math:`\ntau local`",
        default=numpy.array([0.005]),
        domain=Range(lo=0., hi=1000, step=1.0),
        doc="""local tau""")

    cvar = numpy.array([0, 1], dtype=int)

    
    def dfun(self, state_variables, coupling, local_coupling=0.00): 
                
        #number of neurons
        N_e = self.N_tot * (1-self.g)
        N_i = self.N_tot * self.g
        
        #state variables
        
        E, I, C_ee, C_ei, C_ii, W_e, W_i, noise = state_variables
        #delayed state variables
        
        delayed_state = self.history.get_delayed_state(self.local_delay) #50 ms
        E_tau = delayed_state[0, :]
        I_tau = delayed_state[1, :]
   
        
        derivative = numpy.empty_like(state_variables)
        
        # long-range coupling
        c_0 = coupling[0, :]
        
        # short-range (local) coupling
        lc_E = local_coupling * E
        lc_I = local_coupling * I
          
        # external firing rate for the different population
        E_input_excitatory = c_0+lc_E+self.external_input_ex_ex + self.weight_noise * noise
        index_bad_input = numpy.where( E_input_excitatory < 0)
        E_input_excitatory[index_bad_input] = 0.0
        E_input_inhibitory = self.inh_factor * c_0+lc_E+self.external_input_in_ex + self.weight_noise * noise
        index_bad_input = numpy.where( E_input_inhibitory < 0)
        E_input_inhibitory[index_bad_input] = 0.0
        I_input_excitatory = lc_I+self.external_input_ex_in
        I_input_inhibitory = lc_I+self.external_input_in_in
        
        # Transfer function of excitatory and inhibitory neurons with delay
        _TF_e = self.TF_excitatory(E_tau, I_tau, E_input_excitatory, I_input_excitatory, W_e)
        _TF_i = self.TF_inhibitory(E_tau, I_tau, E_input_inhibitory, I_input_inhibitory, W_i)
        
        # Derivatives taken numerically : use a central difference formula with spacing `dx`
        df = 1e-7
        def _diff_fe(TF, fe, fi, fe_ext, fi_ext, W, df=df):
            return (TF(fe+df, fi, fe_ext, fi_ext, W)-TF(fe-df, fi, fe_ext, fi_ext, W))/(2*df*1e3)
        
        def _diff_fi(TF, fe, fi, fe_ext, fi_ext, W, df=df):
            return (TF(fe, fi+df, fe_ext, fi_ext, W)-TF(fe, fi-df, fe_ext, fi_ext, W))/(2*df*1e3)
        
        def _diff2_fe_fe_e(fe, fi, fe_ext, fi_ext, W, df=df):
            TF = self.TF_excitatory
            return (TF(fe+df, fi, fe_ext, fi_ext, W)-2*_TF_e+TF(fe-df, fi, fe_ext, fi_ext, W))/((df*1e3)**2)
        
        def _diff2_fe_fe_i(fe, fi, fe_ext, fi_ext, W, df=df):
            TF = self.TF_inhibitory
            return (TF(fe+df, fi, fe_ext, fi_ext, W)-2*_TF_i+TF(fe-df, fi, fe_ext, fi_ext, W))/((df*1e3)**2)
        
        def _diff2_fi_fe(TF, fe, fi, fe_ext, fi_ext, W, df=df):
            return (_diff_fi(TF, fe+df, fi, fe_ext, fi_ext, W)-_diff_fi(TF, fe-df, fi, fe_ext, fi_ext, W))/(2*df*1e3)
        
        def _diff2_fe_fi(TF, fe, fi, fe_ext, fi_ext, W, df=df):
            return (_diff_fe(TF, fe, fi+df, fe_ext, fi_ext, W)-_diff_fe(TF, fe, fi-df, fe_ext, fi_ext, W))/(2*df*1e3)
        
        def _diff2_fi_fi_e(fe, fi, fe_ext, fi_ext, W, df=df):
            TF = self.TF_excitatory
            return (TF(fe, fi+df, fe_ext, fi_ext, W)-2*_TF_e+TF(fe, fi-df, fe_ext, fi_ext, W))/((df*1e3)**2)
        
        def _diff2_fi_fi_i(fe, fi, fe_ext, fi_ext, W, df=df):
            TF = self.TF_inhibitory
            return (TF(fe, fi+df, fe_ext, fi_ext, W)-2*_TF_i+TF(fe, fi-df, fe_ext, fi_ext, W))/((df*1e3)**2)
        
        _diff_fe_TF_e = _diff_fe(self.TF_excitatory, E_tau, I_tau, E_input_excitatory, I_input_excitatory, W_e)
        _diff_fe_TF_i = _diff_fe(self.TF_inhibitory, E_tau, I_tau, E_input_inhibitory, I_input_inhibitory, W_i)
        _diff_fi_TF_e = _diff_fi(self.TF_excitatory, E_tau, I_tau, E_input_excitatory, I_input_excitatory, W_e)
        _diff_fi_TF_i = _diff_fi(self.TF_inhibitory, E_tau, I_tau, E_input_inhibitory, I_input_inhibitory, W_i)
                
        derivative[0] = (_TF_e - E
                     + .5*C_ee*_diff2_fe_fe_e( E_tau, I_tau, E_input_excitatory, I_input_excitatory, W_e)
                     + .5*C_ei*_diff2_fe_fi(self.TF_excitatory, E_tau, I_tau, E_input_excitatory, I_input_excitatory, W_e)
                     + .5*C_ei*_diff2_fi_fe(self.TF_excitatory, E_tau, I_tau, E_input_excitatory, I_input_excitatory, W_e)
                     + .5*C_ii*_diff2_fi_fi_e( E_tau, I_tau, E_input_excitatory, I_input_excitatory, W_e)
                     )/self.T

        derivative[1] = (_TF_i - I
                     + .5*C_ee*_diff2_fe_fe_i( E_tau, I_tau, E_input_inhibitory, I_input_inhibitory, W_i)
                     + .5*C_ei*_diff2_fe_fi(self.TF_inhibitory, E_tau, I_tau, E_input_inhibitory, I_input_inhibitory, W_i)
                     + .5*C_ei*_diff2_fi_fe(self.TF_inhibitory, E_tau, I_tau, E_input_inhibitory, I_input_inhibitory, W_i)
                     + .5*C_ii*_diff2_fi_fi_i( E_tau, I_tau, E_input_inhibitory, I_input_inhibitory, W_i)
                     )/self.T

        derivative[2] = (_TF_e*(1./self.T-_TF_e)/N_e
                     + (_TF_e-E)**2
                     + 2.*C_ee*_diff_fe_TF_e
                     + 2.*C_ei*_diff_fi_TF_i
                     - 2.*C_ee
                     )/self.T

        derivative[3] = ((_TF_e-E)*(_TF_i-I)
                     + C_ee*_diff_fe_TF_e
                     + C_ei*_diff_fe_TF_i
                     + C_ei*_diff_fi_TF_e
                     + C_ii*_diff_fi_TF_i
                     - 2.*C_ei
                     )/self.T

        derivative[4] = (_TF_i*(1./self.T-_TF_i)/N_i
                     + (_TF_i-I)**2
                     + 2.*C_ii*_diff_fi_TF_i
                     + 2.*C_ei*_diff_fe_TF_e
                     - 2.*C_ii
                     )/self.T
        
        mu_V, sigma_V, T_V = self.get_fluct_regime_vars(
            E_tau, I_tau,
            E_input_excitatory,
            I_input_excitatory,
            W_e, self.Q_e, self.tau_e_e, self.E_e,
            self.Q_i, self.tau_i, self.E_i,
            self.g_L, self.C_m, self.E_L_e, self.N_tot,
            self.p_connect_e, self.p_connect_i, self.g,self.K_ext_e,self.K_ext_i)
        derivative[5] = -W_e/self.tau_w_e+self.b_e*E+self.a_e*(mu_V-self.E_L_e)/self.tau_w_e
        
        mu_V, sigma_V, T_V = self.get_fluct_regime_vars(
            E_tau, I_tau,
            E_input_inhibitory,
            I_input_inhibitory,
            W_i, self.Q_e, self.tau_e_i, self.E_e,
            self.Q_i, self.tau_i, self.E_i,
            self.g_L, self.C_m, self.E_L_i, self.N_tot,
            self.p_connect_e, self.p_connect_i, self.g,self.K_ext_e,self.K_ext_i)
        derivative[6] = -W_i/self.tau_w_i+self.b_i*I+self.a_i*(mu_V-self.E_L_i)/self.tau_w_i
        
        derivative[7] = -noise/self.tau_OU
        
       
        return derivative
    
