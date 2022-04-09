#This file is part of the FOAMS distribution V 0.1 5/7/2022


#import all necessary packages and libraries

import numpy as np
import shelve
import dill
import pickle
import scipy.optimize as optimize


#this is the main function that will be called in the python notebook

def make_metal_fits(elements_to_use, measured_abundances,measured_error,CC_models_to_use, Ia_models_to_use, AGB_models_to_use):
   

    #redfine abundance_func so we can use it here
    
    def abundance_func(x,a,b):
    
        sum_model = a*(CC_model + wind_model + b*Ia_model)

        return np.log10(sum_model[x.astype(int)])
    
    
    
   


    # This loads in all of the variables from the pickle directory made in notebook 2
    
    with open('imf_data_pick.pkl', 'rb') as f:
        data = pickle.load(f)
    
 


    
   

    # This is combining the wind and supernova models of the W18 and the N20 models.
    
    W18_SN_wind_total_yield_elements = data['W18_SN_total_yield_elements'] + data['W18_wind_total_yield_elements']
    data.update({'W18_SN_wind_total_yield_elements': W18_SN_wind_total_yield_elements})
    
    
    N20_SN_wind_total_yield_elements = data['N20_SN_total_yield_elements'] + data['N20_wind_total_yield_elements']
    data.update({'N20_SN_wind_total_yield_elements': N20_SN_wind_total_yield_elements})


    
    # For the Heger 2010 papers, 40 models are available with 10 energies and 4 mixings.
    # Energies available are (index 0 to 9) 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.4, 3.0, 5.0, and 10.0 * 10^51 ergs. 
    # Mixings are 0.0 (no mixing, index 0), 0.001 (standard mixing, index 1), 0.00158 (index 2), and 0.00251 (index 3)

    heger_energy=6   # index 6 = 2.4*10^51 erg
    heger_mixing=0   # index 0 = no mixing

    heger_2010 = data['heger_2010_total_yield_elements'][:,heger_energy,heger_mixing]
    data.update({'heger_2010':heger_2010})

    # Additional models from Heger 2010 may be used by copying the lines above and adding a new model name.
    

    
    all_CC_models = data[CC_models_to_use[0]]


    #This saves all of the models that the user inputs as all_cc_models
    for model in CC_models_to_use[1:]:
        all_CC_models = np.vstack([all_CC_models,data[model]])

    
    #This saves all of the models that the user inputs as all_Ia_models

    all_Ia_models = data[Ia_models_to_use[0]]


    for model in Ia_models_to_use[1:]:
        all_Ia_models = np.vstack([all_Ia_models,data[model]])
    
    
    no_agb = data['karakas_a5_0001_total_yield_elements']*0.
    data.update({'no_agb':no_agb})
    
    
   
    #This saves all of the models that the user inputs as all_wind_models

    all_wind_models = data[AGB_models_to_use[0]]
    


    for model in AGB_models_to_use[1:]:
        all_wind_models = np.vstack([all_wind_models,data[model]])
    
    
    
    #defines normalization, Ia_ratio, and the errors (pcov0 and pcov1) as an array of zeros in the shape of all of the models the user wants to use

    normalization = np.zeros([all_CC_models.shape[0],all_wind_models.shape[0],all_Ia_models.shape[0]])
    Ia_ratio = np.zeros([all_CC_models.shape[0],all_wind_models.shape[0],all_Ia_models.shape[0]])
    pcov_0 = np.zeros([all_CC_models.shape[0],all_wind_models.shape[0],all_Ia_models.shape[0]])
    pcov_1 = np.zeros([all_CC_models.shape[0],all_wind_models.shape[0],all_Ia_models.shape[0]])


    
    #this is an initial guess for the curve_fit function below
    p0 = np.array([1.,0.01])


    # this for loop puts each model the user wants to use into their own arrays to be used in the abundance_func in the curve_fit
    for CC in np.arange(all_CC_models.shape[0]):
        for wind in np.arange(all_wind_models.shape[0]):
            for Ia in np.arange(all_Ia_models.shape[0]):
                CC_model = all_CC_models[CC]
                wind_model = all_wind_models[wind]
                Ia_model = all_Ia_models[Ia]

             

                # Fit solar data with this set of models
                
                #fit to the abundance_func
                try:
                    popt,pcov = optimize.curve_fit(abundance_func,elements_to_use,measured_abundances,p0=p0,sigma=measured_error,absolute_sigma=True,bounds=[[0.,0.],[1.e12,1.]])
                except RuntimeError:
                    print("Error",CC,wind,Ia)

               
                normalization[CC,wind,Ia] = popt[0]
                Ia_ratio[CC,wind,Ia] = popt[1]
                pcov_0[CC,wind,Ia] = pcov[0,0]
                pcov_1[CC,wind,Ia] = pcov[1,1]
             
    #this creates the Ia to Core collapse ratio and the error
    Ia_num = data['Ia_num']
    CC_num = data['cc_num']
                
    Ia_CC_ratio = Ia_ratio*Ia_num/CC_num
    
    Ia_CC_ratio_error = pcov_1 * Ia_CC_ratio/Ia_ratio
    
    print('done')
    
   
    #outputs
    return(normalization,Ia_ratio,Ia_CC_ratio,pcov_0,pcov_1,Ia_CC_ratio_error)

