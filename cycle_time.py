import numpy as np

'''
created by B. N.
2021/01/08

batch_size, arrival_rate, scv_ca, pro_time, scv_pro, setup_time, scv_setup 

'''

class batch_cal ():
    def __init__(self,batch_size, arrival_rate, scv_ca, pro_time, scv_pro, setup_time, scv_setup ):
        self.batch_size = batch_size
        self.arrival_rate = arrival_rate
        self.scv_ca = scv_ca
        self.pro_time = pro_time
        self.scv_pro = scv_pro
        self.setup_time = setup_time
        self.scv_setup = scv_setup

    def withsetup(self):
        eff_pro_time = self.pro_time + self.setup_time / self.batch_size
        ulti = self.arrival_rate * eff_pro_time

        CT = 0

        if ulti < 1:

            scv_cab = self.scv_ca / self.batch_size

                        
            s_sigma_setup = self.scv_setup * self.setup_time**2 #s_sigma_setup = sigma_setup**2

            s_sigma_pro = self.scv_pro * self.pro_time**2  # s_sigma_pro = sigma_pro**2 

            s_eff_sigma = s_sigma_pro + s_sigma_setup/self.batch_size + (self.batch_size - 1)/ self.batch_size**2 * self.setup_time**2

            #s_eff_sigma = eff_sigma**2
            scv_eff_pro = s_eff_sigma / eff_pro_time**2

            scv_eff_pro_b = scv_eff_pro / self.batch_size

            CT_q = (scv_cab + scv_eff_pro_b)/2 * ulti/(1-ulti) * self.batch_size * eff_pro_time

            WTTB = (self.batch_size - 1)/2*(1/self.arrival_rate)

            WTIB = (self.batch_size - 1)/2*self.pro_time

            CT = WTTB + WTIB + CT_q + eff_pro_time

        return CT

if __name__ == "__main__": 
    CT_dict = dict()
    for k in range(3,9):
        func = batch_cal(k, 4, 1.5, 0.15, 0.75, 0.25, 1.5)
        result = func.withsetup()
        CT_dict[k] = result
    print(CT_dict)
    print(min(CT_dict, key=CT_dict.get))