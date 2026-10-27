from parameter_M_Berlin_ARTICLE import Parameter as Parameter_article

class Parameter(Parameter_article):
    def __init__(self, etat='Low'):
        super().__init__(etat=etat)
        # transfer function of Hugo : P_e fitted by Hugo, P_i of Sacha
        self.parameter_model['P_e'] = [-0.04829391053228397, 0.00015312526115162162, -0.0018830942349703766, 0.0022357601238549274, 0.001287915223600697, -0.011256391905113705, 0.007827691831503352, 0.0010805275269078643, -0.003810751044498102, 0.009008257129081343]
        self.parameter_model['P_i'] = [-0.05184978, 0.0061593, -0.01403522, 0.00166511, -0.0020559, 0.00318432, -0.03112775, 0.00656668, 0.00171829, -0.04516385]
