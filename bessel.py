import scipy.special as spe


###################################################################################
###   Script defining spherical Bessel's function and associated Riccati ones   ###
###################################################################################

### hankel's functions 
def hn1(n,x):
    """
    Computation of spherical Bessel's function of the third kind or spherical Hankel's function
    using scipy spherical Bessel's function of the first and second kind.

    Params : 
    --------
                n   : integer - order of the function
                x   : real / array   - param(s) of the function 
    Returns :
    ---------
                output : evaluation of the function in x point(s)

    """
    return spe.spherical_jn(n,x) + 1j*spe.spherical_yn(n,x)


def d_hn1(n,x):
    """
    Computation of spherical Bessel's function of the third kind derivative or spherical Hankel's function derivative
    using scipy spherical Bessel's function of the first and second kind.

    Params : 
    --------
                n   : integer - order of the function
                x   : real / array   - param(s) of the function 
    Returns :
    ---------
                output : evaluation of the function in x point(s)

    """
    return spe.spherical_jn(n,x,derivative =True) + 1j*spe.spherical_yn(n,x,derivative =True)



### riccati bessel functions

def rj(n,x):

    return x* spe.spherical_jn(n,x)


def d_rj(n,x):

    return spe.spherical_jn(n,x) + x*spe.spherical_jn(n,x, derivative=True)

def rh(n,x):

    return x*hn1(n,x)

def d_rh(n,x):

    return hn1(n,x) + x*d_hn1(n,x)




if __name__=="__main__":
    

    print("->    little test !!!!\n")
    print("j01(3)")
    print(spe.spherical_jn(0,3))
    print("---------------\n")


    print("->    little test !!!!\n")
    print("h01(3)")
    print(hn1(0,3))
    print("---------------\n")

    print("d_h01(3)")
    print(d_hn1(0,3))
    print("---------------\n")

    print("rh_01(3)")
    print(rh(0,3))
    print("VS \n riccati bessel of python :")
    print(spe.riccati_jn(0,3)[0] + 1j*spe.riccati_yn(0,3)[0])
    print("---------------\n")

    print("d_rh_01(3)")
    print(d_rh(0,3))
    print("VS \n riccati bessel of python :")
    print(spe.riccati_jn(0,3)[1] + 1j*spe.riccati_yn(0,3)[1])



    print("\nrj_01(3)")
    print(rj(0,3))
    print("VS \n riccati bessel of python :")
    print(spe.riccati_jn(0,3)[0] )

    print("\nd_rj_01(3)")
    print(d_rj(0,3))
    print("VS \n riccati bessel of python :")
    print(spe.riccati_jn(0,3)[1] )




    print("---------------\n")




