
# import numpy
import numpy as np


# define a function for the powerlaw imf
def make_imf_powerlaw(n, m_min, m_max, m_cut, imf_index ):
    
    index = imf_index
    a = np.arange(n+1)/n

    b = 10.**((np.log10(m_max)-np.log10(m_min))*a) * m_min

    mass_limits = b
    
    c = b**index
    b0 = b[0:n]
    b1 = b[1:(n+1)]
    
    d = (1./(index+1) * b1**(index+1) - 1./(index+1) * b0**(index+1))
    num_m_1 = d
        
    
    e = (1./(index+2) * b1**(index+2) - 1./(index+2) * b0**(index+2))
    f = e/num_m_1
    
    mass = f
        
    z0 = np.zeros(n)
    
    num_m_2 = np.where(mass > m_cut,z0,num_m_1)
    
    total_mass = np.sum(num_m_2*mass)
    num_m = num_m_2/total_mass
    
    #print(num_m)
    #print(mass)
    #print(np.sum(num_m*mass))
    
  
    return mass, num_m, mass_limits


#define a function for the kroupa imf
def make_imf_kroupa(n, m_min, m_max, m_cut, imf_index):
    index = imf_index

    index1 = -0.3
    index2 = -1.3
    
    mass1 = 0.08
    mass2 = 0.5
    
    a = np.arange(n+1)/n
    
    b = 10.**((np.log10(m_max)-np.log10(m_min))*a) * m_min
        
    mass_limits = b
        
    c = np.zeros(n+1)
    
    c = np.where(b < mass1,b**index1,c)
    
    c = np.where((b >= mass1) & (b < mass2) , b**index2 * mass1**(index1)/(mass1**index2),c)
        
    c = np.where(b >= mass2, b**index * mass1**index1/mass1**index2 * mass2**index2/mass2**index,c)
    
    b0 = b[0:n]
    b1 = b[1:(n+1)]
    
    d = (1./(index+1) * b1**(index+1) - 1./(index+1) * b0**(index+1))
    num_m_1 = d
    
    d = c*b
    num_m_1 = d[0:n]
    
    e = (1./(index+2) * b1**(index+2) - 1./(index+2) * b0**(index+2))
    f = e/num_m_1
        
    mass = f
    
    mass = b[0:n]
    
    
    z0 = np.zeros(n)
    
    num_m_2 = np.where(mass > m_cut,z0,num_m_1)
    
    total_mass = np.sum(num_m_2*mass)
    num_m = num_m_2/total_mass
    
    return mass, num_m, mass_limits




#define a function for the yield per bin and star for each model
def model_yield(type2_min_mass,type2_max_mass, mass, num_m, star_masses, isotope_names, big_data):
    
    n = np.size(mass)
    
#    yield_per_star = np.zeros((n,len(isotope_names)-1))
    yield_per_star = np.zeros((n,len(isotope_names)))

    print(yield_per_star.shape)
    print(big_data.shape)

    for i in np.arange(n):
        if mass[i] < type2_min_mass:
            yield_per_star[i,:] = 0.
        if ((mass[i] >= type2_min_mass) & (mass[i]<=type2_max_mass)):
            #figure out which SN data to use
            # find index of clostest mass of W18 models.  
            mass_of_star = mass[i]
            index_star = (np.abs(star_masses - mass_of_star)).argmin()
        #index_star = 1
            print('index_star = ',index_star,star_masses[index_star])
            yield_per_star[i,:] = big_data[index_star,:] #W18_yield[index_star,:]
    
        if mass[i] > type2_max_mass:
            yield_per_star[i,:] = 0.

    print('yield_per_star[:,13] = ',yield_per_star[:,13])
    print('yield_per_star[:,61] = ',yield_per_star[:,61])
    print('yield_per_star O/Fe = ',yield_per_star[:,13]/yield_per_star[:,61])
    
    
    
    yield_per_bin = yield_per_star

    for i in np.arange(n):
        yield_per_bin[i,:] = num_m[i]*yield_per_star[i,:]

    return yield_per_bin, yield_per_star


