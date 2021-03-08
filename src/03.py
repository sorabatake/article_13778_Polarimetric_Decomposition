import os
import struct
import numpy as np
import glob
from tqdm import tqdm
import pickle
import cv2
import matplotlib.pyplot as plt
from slcinfo import SLC_L11
import gc

# Fields
SAVE_DIRECTORY="./data/"
def slc2cov(_hh, _hv, _vv, _vh, _nb=2, _cross=0):    
    [m,n] = _hh.shape
    result = np.zeros((m,n,9),dtype=np.complex64)
    if _cross == 0:
        _cross = (_hv + _vh) / 2
    else:
        _cross = _hv   
    cross = cross * np.sqrt(2)
    for i in tqdm(np.arange(0, m)):
        for j in np.arange(0, n):
            
            i_one = np.array((0, i - _nb))
            i_two = np.array((i+_nb+1, m))
            
            j_one = np.array((0, j-_nb))
            j_two = np.array((j+_nb + 1, n))
    
            hh_tmp = np.copy(_hh[np.max(i_one):np.min(i_two),np.max(j_one):np.min(j_two)])
            vv_tmp = np.copy(_vv[np.max(i_one):np.min(i_two),np.max(j_one):np.min(j_two)])
            cross_tmp = np.copy(cross[np.max(i_one):np.min(i_two),np.max(j_one):np.min(j_two)])

            [mbis,nbis]=hh_tmp.shape
            
            hh_tmp = hh_tmp.reshape((mbis*nbis, 1))
            vv_tmp = vv_tmp.reshape((mbis*nbis, 1))
            cross_tmp = cross_tmp.reshape((mbis*nbis, 1))

            data = np.hstack((hh_tmp, cross_tmp, vv_tmp))
            
            cov = np.dot(data.T,np.conjugate(data))/(mbis*nbis)
            result[i,j,:] = cov.reshape(9,)   
    return result

def cov2fdd(data):
    [m, n, p] = data.shape
    pv_array = np.zeros((m, n))
    pd_array = np.zeros((m, n))
    ps_array = np.zeros((m, n))    
    for i in tqdm(np.arange(0, m)):
        for j in np.arange(0, n):        
            cov=(np.copy(data[i,j,:])).reshape((3, 3))
            fv=(3*cov[1,1]/2).real            
            if (cov[0,2]).real<=0:
                beta=1
                
                Af=(cov[0,0]-fv).real
                Bf=(cov[2,2]-fv).real
                Cf=((cov[0,2]).real-fv/3).real
                Df=(cov[0,2]).imag
                
                if (Af+Bf-2*Cf)!=0:
                    fd=(Df**2+(Cf-Bf)**2)/(Af+Bf-2*Cf)
                else:
                    fd=Bf
                
                fs=Bf-fd
                
                if fd!=0:
                    alpha=np.complex((Cf-Bf)/fd+1, Df/fd)
                else:
                    alpha=np.complex(1,1)
                
            else:           
                alpha=-1
                
                Af=(cov[0,0]-fv).real
                Bf=(cov[2,2]-fv).real
                Cf=(cov[0,2]).real-fv/3
                
                if (Af+Bf+2*Cf)!=0:                    
                    fs=(Cf+Bf)**2/(Af+Bf+2*Cf)
                else:
                    fs=Bf
                    
                fd=Bf-fs
                beta=(Af+Cf)/(Cf+Bf)       
                
            pv=(8*fv/3)
            ps=(fs*(1+beta*np.conjugate(beta))).real
            pd=(fd*(1+alpha*np.conjugate(alpha))).real
            
            ps_array[i, j] = ps
            pd_array[i, j] = pd
            pv_array[i, j] = pv
    return np.asarray([ps_array, pd_array, pv_array])

def get_intensity(_slc): 
    return 20*np.log10(abs(_slc)) - 83.0 - 32.0

def normalization(_x):
    return (_x-np.amin(_x))/(np.amax(_x)-np.amin(_x))

def adjust_value(_x):
    _x = np.asarray(normalization(_x)*255, dtype="uint8")
    return cv2.equalizeHist(_x)

def slc2pauli(_hh, _hv, _vv, _vh, _cross=0):
    if _cross == 0:
        _cross = (_hv + _vh) / 2
    else:
        _cross = _hv    
    single = ((_hh+_vv)/np.sqrt(2)) # single, odd-bounce.
    double = ((_hh-_vv)/np.sqrt(2)) # double, even-bounce
    vol = np.sqrt(2)*_cross # vol
    return np.asarray([single, double, vol])

def grayplot(_file_name, _v):
    plt.figure()
    plt.imshow(_v, cmap = "gray")
    plt.imsave(_file_name, _v, cmap = "gray")
    print("[Save]", _file_name)
    
def rgbplot(_file_name, _alpha_r, _gamma_g, _beta_b):
    plttot=(np.dstack([_alpha_r, _gamma_g,_beta_b]))
    plttot=np.asarray(plttot, dtype="uint8")
    plt.figure()
    plt.imshow(plttot)
    plt.imsave(_file_name, plttot)
    print("[Save]", _file_name)
    return

# Entry point
def main():
    offset_x = 3000
    offset_y = 8000
    limit_length = 2000
    slc_hh = 0
    slc_hv = 0
    slc_vv = 0
    slc_vh = 0
    product_list = os.listdir(SAVE_DIRECTORY)
    for product_name in product_list:
        slc_list = glob.glob(SAVE_DIRECTORY + product_name + "/IMG-*.pkl")
        for slc_name in slc_list:
            with open(slc_name, "rb") as f:
                slc = pickle.load(f)
            if "IMG-HV" in slc_name:
                slc_hv = slc.data
            if "IMG-HH" in slc_name:
                slc_hh = slc.data
            if "IMG-VV" in slc_name:
                slc_vv = slc.data
            if "IMG-VH" in slc_name:
                slc_vh = slc.data
            slc_tmp = slc2intensity(slc.data)
            grayplot(slc_name + "_p.png", adjust_value(slc_tmp))
            slc = None
            gc.collect()
        
        ##################################
        # pauli
        ##################################
        result = slc2pauli(_hh=slc_hh, _hv=slc_hv, _vv=slc_vv, _vh=slc_vh) # alpha, beta, gmma
        [s, d, v] = get_intensity(result)
        rgbplot(SAVE_DIRECTORY + product_name + "/Pauli_decomposition.png", _alpha_r=adjust_value(d), _gamma_g=adjust_value(v), _beta_b=adjust_value(s))
        ##################################
        # freeman-duran
        ##################################
        # convert slc to cov
        cov_file_name = SAVE_DIRECTORY + product_name + "/cov.pkl" 
        if not os.path.exists(cov_file_name):
            cov = slc2cov(slc_hh, slc_vv, slc_vv, slc_vh, 2)
            with open(cov_file_name, "wb") as f:
                pickle.dump(cov, f)
        with open(cov_file_name, "rb") as f:
            cov = pickle.load(f)

        # convert freeman-duran
        freeman_file_name = SAVE_DIRECTORY + product_name + "/fdd.pkl" 
        if not os.path.exists(freeman_file_name):
            fdd = cov2fdd(cov)
            with open(freeman_file_name, "wb") as f:
                pickle.dump(fdd, f)
        with open(freeman_file_name, "rb") as f:
            fdd = pickle.load(f)
        [s, d, v] = get_intensity(fdd)
        rgbplot(SAVE_DIRECTORY + product_name + "/Freeman_decomposition.png", _alpha_r=adjust_value(d), _gamma_g=adjust_value(v), _beta_b=adjust_value(s))

if __name__=="__main__":
       main()

