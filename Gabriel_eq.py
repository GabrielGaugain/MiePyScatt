import numpy as np 
import pandas as pd
import ntpath
import matplotlib.pyplot as plt
from scipy.constants import mu_0, epsilon_0

####################################################################
###   Script for dielectric properties of tissues calculations   ###
####################################################################



def Gabriels_Permitivity(omega, ef, sigma, deltas, alphas, taus):
    """
    Implementation of the Gabriel's equation from S Gabriel and R W Lau and C Gabriel, 1996 
    "The dielectric properties of biological tissues:  III. Parametric models for the dielectric spectrum of tissues "
    to compute permitivity of tissues.
    
    Params : 
    --------
                omega : frequency | real 
                ef    : permtivity at HF | real
                sigma : conductivity | real
                deltas,alphas , taus : lists of 4 parameters in the sum of Gabriel's equation (Cole Cole interaction) |List[real]

    Returns :
    ---------

                espilon :  permitivity | complex 

    """

    ## Only if there are 4 parameters in each list for Alphas, taus and deltas
    assert(deltas.size == 4)
    assert(alphas.size == 4)
    assert(taus.size == 4)

    ## Declaration of the permitivity and init at ef + cond contribution 
    epsilon = ef + sigma/(1j*omega*epsilon_0 )

    ## loop over the list of params = sum in the equation
    for i in range(4):
        epsilon = epsilon + deltas[i]/( 1 + ( 1j*omega*taus[i])**(1-alphas[i]) )

    ## return the resulting complex electic permitivity
    return epsilon

## ---- END Gabriels_Permitivity() ---- ##


def diel_read(path="", filename ="diel.xlsx"):
    """
    Open in read mode the excel file containing all parameters needed to compute the permitivity.
    Using pandas library to read the sheet.

    Params :
    --------
            path :  the path where the file is located. Default = current path (ie same path as the script)
            filename :  name of the excel file containing params. Default = "diel.xlsx"

    Returns :
    ---------
            The pandas structure (array containing info of the file)

    """

    try : 
        # reading the sheet with pandas function read_excel
        if(path != ""):
            dataframe = pd.read_excel( ntpath.join(path,filename) )
        else :
            dataframe = pd.read_excel(filename)
        
        return dataframe

    except:
        print("ERROR : Could not read the file (not fount or cant read it) " )

## ---- END diel_read() ---- ##


def get_Params(df, tissues):
    
    if isinstance(tissues, str):
        tissues = np.array([tissues])

    try:
        N_tissues = np.size(tissues)

        Params = {  "tissues" : tissues, \
                    "ef":np.zeros((N_tissues,1)), "sig":np.zeros((N_tissues,1)), \
                    "alphas":np.zeros((N_tissues,4)), "taus": np.zeros((N_tissues,4)), "deltas": np.zeros((N_tissues,4)) \
                }


        for i,tissue in enumerate(tissues):

            ## selecting the tissue of interest to compute its permitivity
            TOI = df[df["Name"]== tissue]
            # print(TOI)  # Uncomment to print the line corresponding to the tissue of interest

            ## extracting the params of this tissues from the dataframe
            alphas = TOI.iloc[0, [ ('alf' in col) for col in TOI.columns ]].values
            Params["alphas"][i,:] = alphas

            taus = TOI.iloc[0, [('tau' in col) for col in TOI.columns] ].values * np.array([ 10**(-12), 10**(-9), 10**(-6), 10**(-3)] )    ## Careful to unity for time variables (ps, ns, µs, ms)
            Params["taus"][i,:] = taus

            deltas = TOI.iloc[0, [('del' in col) for col in TOI.columns] ].values
            Params["deltas"][i,:] = deltas

            ef = TOI["ef"].values[0]
            Params["ef"][i] = ef

            sigma = TOI["sig"].values[0]
            Params["sig"][i] = sigma
        
    except:
        print("WARNING : the tissue is not find in the table, use a tissue contained in the file")
        Params = []

    return Params

## ---- END get_Params() ---- ##

def diel_save(dict):
    """
    Function to save Permitivity data in an excel file 

    Params :
    --------

    Returns :
    ---------

            isSaved :   Boolean | True if the data is correctly saved, else False

    """


    ## SECTION TO SAVE THE DATA in an excel file
    """
    ## create pandas dataframes with first col as names and others perm at different freq in Hz
    Real_perm = pd.DataFrame(np.concatenate(( names, np.real(perm)), axis=1  ).T, 
                            
                           )

    #print( np.concatenate(( names, np.real(perm)[0]), axis=1  ) )
  
    Imag_perm = pd.DataFrame(np.concatenate(( names , np.imag(perm)), axis=1  ), 
                             columns=['Name']+freqs
                            )    
    

    ## Write the data in differents sheets for real and imaginary parts
    
    writer = pd.ExcelWriter('PermitivitiesTest.xlsx')

    Real_perm.to_excel(writer, 'Real_eps')
    #Imag_perm.to_excel(writer, 'Imag eps')

    writer.save()
    """    

    ## To be completed
    isSaved = False

    return isSaved

## ---- END diel_save() ---- ##



def spectrum_computing(Params, f= np.logspace(1,8,10**4) , tissue = "some tissue", plot=False ):
    """
    Compute the Gabriel's equation for the permitivity between 10hz to 100 kHz (default) or the frequency range provided

    Params :
    --------

                params  :   dictionnary of the parameters needed to compute the permitivity. looks like : {"ef": 0, "sig": 0 ,"alphas":np.array([]), "taus": np.array([]), "deltas": np.array([]) }
                f       :   array of the frequecies on which to compute the permitivity |   DEFAULT = [ 10Hz, 100kHz] pas de 1 Hz
                tissue  :   name of the tissue investigated (just for the title of the plot) |   DEFAULT = "some tissue"
                plot    :   Boolean, true to plot the real part of perm and false to just return the array

    Returns :
    ---------
                permitivities   :   complex array containing permitivity in the range 10Hz-100kHz

    """   

    ## init array
    permitivities = np.zeros(np.size(f))

    ## Compute the perm thanks to the function Gabriels_Permitivity. 
    permitivities = Gabriels_Permitivity( 2*np.pi*f, Params["ef"], Params["sig"], Params["deltas"], Params["alphas"], Params["taus"])

    ## Get conductivitie : \omega Im(e)e_0
    conduct = -np.imag(permitivities)*epsilon_0*2*np.pi*f

    if plot:
        ## plot of the real part of relative permitivity in log-log scale
        fig, ax1 = plt.subplots()
        plt.title("Dielectric properties of  "+tissue)
        ax2 = ax1.twinx() ## to add another y-axis for conductivity

        color1, color2 = "darkblue", "darkred"

        ax1.loglog(f, permitivities.real, color=color1)
        ax1.set_ylabel("Relative permitivity ", color=color1)
        ax1.tick_params(axis='y', labelcolor=color1)
        ax1.set_xlabel( "frequency (Hz)")

        ax2.loglog(f, conduct, color=color2, linestyle='dashed' )
        ax2.set_ylabel("Conductivity (S/m)", color=color2)
        ax2.tick_params(axis='y', labelcolor=color2)
        
        fig.tight_layout()
        #plt.xlim((5, 10000000))
        #plt.ylim((permitivities[-1].real/10,permitivities[0].real*10 )) ## set to have a little blanc range on top and bottom of the curve
        plt.grid()
        plt.show()


    return permitivities

## ---- END spectrum_computing() ---- ##

def get_allP_to_OneW(f, dataframe ):
    """
    Function to get the permitivity of all tissues presented in the file at one frequency.

    Params :
    --------

                f       :   frequency (Hz) in which to compute permitivities
                dataframe   :   dataframe extracted from the xls file

    Returns : 
    ---------

                eps         :   Permitivities of all (N) tissues at the frequency omega

    """
    
    ## get all tissues (names)
    tissues = dataframe["Name"][:]
    N = tissues.size

    
    if isinstance(f, list) : f = np.array(f).T
    
    if isinstance(f, np.ndarray) : 

         ## init the array of permitivities
        eps = np.empty((N, f.size), dtype=complex)      

    else:    
         ## init the array of permitivities
        eps = np.empty(N, dtype=complex)



    for i in range(N):

        ## selecting the tissue of interest to compute its permitivity
        TOI = dataframe[dataframe["Name"]== tissues[i]]

        ## extracting the params of this tissues from the dataframe
        alphas = TOI.iloc[0, [ ('alf' in col) for col in TOI.columns ]].values
        taus = TOI.iloc[0, [('tau' in col) for col in TOI.columns] ].values * np.array([ 10**(-12), 10**(-9), 10**(-6), 10**(-3)] )    ## Careful to unity for time variables (ps, ns, µs, ms)
        deltas = TOI.iloc[0, [('del' in col) for col in TOI.columns] ].values
        ef = TOI["ef"].values[0]
        sigma = TOI["sig"].values[0]

        ##
        if isinstance(f, np.ndarray) : 
            eps[i,:] = Gabriels_Permitivity( 2*np.pi*f , ef ,sigma , deltas, alphas, taus)
        else:
            eps[i] = Gabriels_Permitivity( 2*np.pi*f , ef ,sigma , deltas, alphas, taus)

    return eps

## ---- END get_allP_to_OneW() ---- ##


################ => most usefull function to get properties for another script

def get_diel_properties(tissues, freqs= np.logspace(1,8,10**5)):
    """
    Compute the Gabriel's equation for the permitivity between 10hz to 10 MHz (default) or the frequency range provided,
    and returns the relative permitivities (real part) and conductivities for each tissue.

    Params :
    --------
                tissues :   str or array[str] | str name of the tissue investigated (just for the title of the plot) |   DEFAULT = "some tissue"
                f       :   array[float]      | frequecies on which to compute the permitivity |   DEFAULT = [ 10Hz, 100kHz] pas de 1 Hz

    Returns :
    ---------
                permitivities   :   array[float]  | Relative permitivities of each tissues computed of the provided freq range
                conductivities  :   array[float]  | Relative permitivities of each tissues computed of the provided freq range

    """   

    if isinstance(tissues, str):
        tissues = np.array([tissues])

    ## read the file
    df = diel_read()

    ## get parameters needed for Gabriel's model (Cole-Cole)S
    P = get_Params(df, tissues)


    if (P == []):
        print("Failed to get parameters from file - the provided tissues do not exist")
        return 

    permitivities, conductivities = np.zeros( (np.size(tissues), np.size(freqs)) ), np.zeros( (np.size(tissues), np.size(freqs)) )
    
    for i,tissue in enumerate(tissues):

        ## get parameters for the i th tissue
        Params = {"ef": P["ef"][i], "sig":P["sig"][i],"alphas":P["alphas"][i,:], "taus": P["taus"][i,:], "deltas": P["deltas"][i,:]} 

        ## computing complex permitivity and plot it (10Hz - 10 kHz by default) 
        spec = spectrum_computing(Params,f=freqs, tissue=tissue, plot=False)

        permitivities[i, :] = spec.real
        conductivities[i,:] =    -np.imag(spec)*epsilon_0*2*np.pi*freqs


    return permitivities , conductivities

## ---- END get_diel_properties() ---- ##



def main():
    """
    Main function of this script => to modify to play with functions and get the result you want
    """

    ## tissue(s) to investigate / Tissue(s) of Interest
    #tissue = "Brain (Grey Matter)"
    tissues = np.array(["Bone (Cortical)"])
    #tissues = np.array(["Bone (Cortical)", "Brain (White Matter)", "Brain (Grey Matter)", "Skin","Cerebrospinal Fluid"])

    ## read the file 
    df = diel_read()
    #print(df.head()) # Uncomment to print the first lines of the dataframe in the console

    ## array to store the frequency(ies) at which the relative error between QS and FW becomes higher than 1%
    min_f = np.zeros(tissues.shape)

    ## get parameters needed for Gabriel's model (Cole-Cole)S
    P = get_Params(df, tissues)


    for i,tissue in enumerate(tissues):

        ## get parameters for the i th tissue
        Params = {"ef": P["ef"][i], "sig":P["sig"][i],"alphas":P["alphas"][i,:], "taus": P["taus"][i,:], "deltas": P["deltas"][i,:]} 
        
        ## computing complex permitivity and plot it (10Hz - 100 kHz by default) 
        freqs= np.logspace(1,8,10**5)
        spec = spectrum_computing(Params,f=freqs, tissue=tissue, plot=True)

        Rel_perm = spec.real
        #print(spec.size)

        ####  Plotting results #################################

        ## plot of the real part of permitivity in log-log scales
        plt.figure()
        plt.subplot(211)
        plt.loglog(freqs, Rel_perm, 'k')
        plt.xlabel( "frequency (Hz)")
        plt.ylabel("Relative permitivity ")
        plt.title("Relative permitivity of "+  tissue)
        #plt.xlim((5, 10000000))
        #plt.ylim((permitivities[-1].real/10,permitivities[0].real*10 )) ## set to have a little blanc range on top and bottom of the curve
        plt.grid()


        ### Inspection of the caracteristique error factor 
        plt.subplot(212)
        #plt.semilogx(freqs, freqs**2 *Rel_perm *0.16**2 )
        #plt.semilogx(freqs, np.ones(Rel_perm.shape)*0.01)
        plt.loglog(freqs, freqs**2 *Rel_perm *0.16**2 *mu_0*epsilon_0 , 'k-')
        plt.loglog(freqs, np.ones(Rel_perm.shape)*0.01, linestyle ="--")
        plt.xlabel( "frequency (Hz)")
        plt.ylabel( r'$ \omega^2  \mu \epsilon $')
        plt.show()  
    
        ## get all freq for which the factor is greater than 1% and then get the first value (limit)
        err = freqs[  freqs**2 *Rel_perm *0.16**2 *mu_0*epsilon_0 > np.ones(Rel_perm.shape)*0.01   ]
        min_f[i] = err[0]

    ## end for tissues

    print(tissues)
    print(min_f)
    print(np.log10(min_f))



    ##  SECTION TO SAVE CONDUCT AND PERM IN TXT FILE FOR COMSOL IMPORT
    """
    with open("RelPerm_"+tissue+".txt", 'w') as f_perm, open("Conduct_"+tissue+".txt", 'w') as f_cond:
        for i in range(spec.size):
            f_perm.write( str(freqs[i]) + ' '  +str(Rel_perm[i])+ '\n')
            f_cond.write( str(freqs[i]) + ' '  + str(conduct[i])+ '\n')
    """

    ## SECTION COMPUTING ALL perm for ALL tissues
    """    
    freqs = [10,50,100,200,500,1000,2500,10000]
    ## freqs = np.concatenate(( np.arange(0,1000,10) , np.arange(1250,100000,250) , np.arange(101000,10000000,1000) ))
    perm  =get_allP_to_OneW(  freqs, df ) 
    print(perm.shape)

    ## getting name of tissues
    names = np.array( "E_" + df["Name"][:])
    names = names.reshape((names.size,1))

    e_max = np.nanmax( np.real(perm), axis=0 )
    """


## ---- END main() ---- ##



if __name__ == "__main__":
    # Call the main function 
    main()





