"""Collection of class definitions and functions for time-temperature superposition."""



class william_landel_ferry_model():
    
    def __init__(self, C1, C2, T_ref:float) -> None:
        
        if not isinstance(C1, (int,float)) or not isinstance(C2, (int,float)) or not isinstance(T_ref, (int,float)):
            raise TypeError('Empircal constants C1 and C2, and reference temperature T_ref, must be a numerical type.')
        
        self.C1 = C1
        self.C2 = C2
        self.T_ref = T_ref
        
        

wlf_model = lambda C_1, C_2, T_ref, T: -C_1 * (T - T_ref) / (C_2 + (T - T_ref))
