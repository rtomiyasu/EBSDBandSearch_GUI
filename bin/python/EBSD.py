# -*- coding: utf-8 -*-
import numpy as np
import time
#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
import sys
import os, os.path
import importlib
import bisect
import pickle
from skimage.io import imread, imsave
from skimage.transform import radon, rescale
from skimage import exposure

import mask

os.chdir(os.path.dirname(os.path.abspath(__file__)))

# (rho, theta)に対応する直線の、画像上の端点を算出．mask.pyのgetLine関数と基本的に同じだが、
# Plotする際は、(x,y) -> (y-image_shape[0]+1, x)の座標変換が必要．
def getLineForDisplay(image_shape, rho, theta, LineX, LineY):
    LineX.clear()
    LineY.clear()
    image_o = [image_shape[0]//2,image_shape[1]//2]
    # 画像のエッジ x = -0.5, x = image.shape[1]-0.5との交点を求める.
    if theta != 0.0:
        x = np.array([-0.5, image_shape[1]-0.5])
        y = image_shape[0]-1-image_o[0] + (x - image_o[1])/np.tan(theta) - rho/np.sin(theta) 
        for i in np.arange(0,2):
            if 0 <= y[i] and y[i] <= image_shape[0]:
                LineX.append(x[i])
                LineY.append(y[i])
    
    # 画像のエッジ y = -0.5, y = image.shape[0]-0.5との交点を求める
    if theta != np.pi*0.5:
        y = np.array([-0.5, image_shape[0]-0.5])
        x = image_o[1] + (y - image_shape[0]+1+image_o[0])*np.tan(theta) + rho/np.cos(theta) 
        for i in np.arange(0,2):
            if 0 <= x[i] and x[i] <= image_shape[1]:
                LineX.append(x[i])
                LineY.append(y[i])

# 配列ArrayRhoThetaに格納された各ArrayRhoTheta[i]に対応する、画像上の端点を算出．
# ただし、各ArrayRhoTheta[i]は(rho, theta)に対応する配列のindex(整数値)とする。
# 上記のgetLine関数を使用．
# def getLines_irho_itheta(theta, image_shape, rho_o, ArrayRhoTheta, ArrayLines):
#     ArrayLines.clear()
#     for rt in ArrayRhoTheta:
#         LineX = []
#         LineY = []
#         getLine(image_shape, rt[0]-rho_o, np.deg2rad(theta[rt[1]]), LineX, LineY)
#         ArrayLines.append([LineX, LineY])


# 多項式フィッティングによる、平滑化・微分値計算に用いる:
def putX(x_,n,num_points):
    return x_[n-num_points:n+num_points+1] - x_[n]

def putY(y_,n,num_points):
    return y_[n-num_points:n+num_points+1]

def putCoef(x_,y_,n,num_points,deg):
    return np.polyfit(putX(x_,n,num_points),putY(y_,n,num_points),deg)

# ラドン変換の誤差値(各(rho, sigma)に対応する直線の画像内線分の長さ*sigma)を求める
def calRadonGosa(thetas, image, sinogram, sigma, circle, ArraySinogramErrors):
    ArraySinogramErrors.clear()
    rho_o = sinogram.shape[0]//2
    for n in range(sinogram.shape[0]): # index for rho
        lengths = []
        for m in range(sinogram.shape[1]): # index for theta 
            # まず各(rho, theta)に対応する線分の長さを求める
            LineX = []
            LineY = []
            # getLines_irho_itheta(theta, image.shape, sinogram.shape[0]//2, [[n,m]], ls)
            mask.getLine(image.shape, n-rho_o, np.deg2rad(thetas[m]), LineX, LineY, circle)
            llength = 0.0
            for i in range(1, len(LineX)):
                # ls[0][0] : x-coordinates, ls[0][1] : y-coordinates, 
                llength = max(llength, (LineX[0]-LineX[i])**2 + (LineY[0]-LineY[i])**2)
            lengths.append( np.sqrt(np.sqrt(llength))*sigma )
        ArraySinogramErrors.append(lengths)

# 球面座標(cos(sigma)*cos(phi),cos(sigma)*sin(phi),sin(sigma))を格納: 
# ただしx、ｙ軸は画像上を右・上向きに、z軸は画像に垂直に取る (この後のindexingで用いている座標系)。   
class SphericalCoordinate:
    def __init__(self):
        self.phi = 0.0   # radian
        self.sigma = 0.0 # radian
    
# 検出したピーク座標を格納: 
class PeakData:
    def __init__(self):
        self.iinterval = [0,0]     # 区間(整数値). A peak exists in [x[iinterval[0]], x[iinterval[1]]).
        self.dinterval = [0.0,0.0] # 区間(実数値). A peak exists in [dinterval[0], dinterval[1]).
        self.peak_height = 0.0     # Hough変換におけるピーク高さ（推定値）
        
# 上向き下向きのピークのペアが片側のバンドエッジに対応する． 
class EdgeData:
    def __init__(self):
        self.left = PeakData()
        self.right = PeakData()
        self.zero_rho = 0.
        self.isPlusMinus = False # 0:-+型, 1:+-型 (rho方向の2次微分に関して).

# バンドエッジ, バンドセンター座標に関する情報を格納
class BandData:
    def __init__(self):
        # self.itheta = 0
        self.edge_rhos = [0.0, 0.0] # edge_rhos[0] <= edge_rhos[1] 
        self.center_rt = [0.0, 0.0] # (rho, theta (deg.)) of the center line.
        self.center_SC = SphericalCoordinate() # spherical coordinate of the center line.
        self.BraggAngle = 0.0 # in radian
        
        # Parameters used for comparing bands.
        self.dark_edge_rhos = [0.0, 0.0] # バンドエッジが同じrho座標にあるかどうかを判定するために幅を持たせておく。
        self.convolution = 0.0
        # self.num_band_pairs_good_match = []
    def putTheta(self):
        return self.center_rt[1]
    def putConvolution(self):
        return self.convolution
        # return self.dark_edge_heights[0] + self.dark_edge_heights[1]
    def putEdgeRanges(self):
        return [[self.dark_edge_rhos[0], self.edge_rhos[0]], [self.edge_rhos[1], self.dark_edge_rhos[1]]]
    def setEdges(self, theta, edge_range1, edge_range2):
        self.center_rt[1] = theta
        self.edge_rhos = [edge_range1[1], edge_range2[0]]
        self.dark_edge_rhos = [edge_range1[0], edge_range2[1]]
    # edge_rhosを用いて center_rt, center_sc, BraggAngleをsetする．
    def setCenter(self, PC, image_o, center_rho=None):
        [self.center_SC, self.BraggAngle] = transform_RhoTheta_to_SphericalCordinate(image_o, PC, self.putTheta(), self.edge_rhos[0], self.edge_rhos[1], center_rho)
        self.center_rt = transform_SphericalCordinate_to_RhoTheta(image_o, PC, self.center_SC)
        self.center_rt[1] = np.rad2deg(self.center_rt[1]) # band center　座標()のthetaをradianをdegreeに変換

# バンドセンターrhoC (rhoC=Noneの場合, エッジ(rho1, theta), (rho2, theta)に関する情報)から、
# バンドセンターにPCから下した垂線の足の球面座標よる表示SCとバンド幅に関する情報を得る。
def transform_RhoTheta_to_SphericalCordinate(image_o, PC, theta, rho1, rho2, rhoC=None):
    SC = SphericalCoordinate()
    SC.phi = np.deg2rad(theta)
    # rhoについては下式:
    # rho = (PC[1] - image_o[1])*cos(SC.phi) - (PC[0] - image_o[0])*sin(SC.phi) + tan(SC.sigma)*PC[2] 
    rho0 = (PC[1] - image_o[1])*np.cos(SC.phi) - (PC[0] - image_o[0])*np.sin(SC.phi)
    sigma1 = np.arctan((rho1 - rho0)/PC[2])
    sigma2 = np.arctan((rho2 - rho0)/PC[2])
    SC.sigma = (sigma1 + sigma2)*0.5 if rhoC==None else np.arctan((rhoC - rho0)/PC[2])
    bangle = abs(sigma1 - sigma2)*0.5
    if SC.sigma < 0.0:
        SC.phi += np.pi
        SC.sigma *= -1.
    return [SC, bangle]

# バンドセンターの球面座標SCによる表示から、(rho,theta)による表示を得る
def transform_SphericalCordinate_to_RhoTheta(image_o, PC, SC):
    # 以下は同じ直線:
    # sin(theta)*(y - image_o[1]) + cos(theta)*(x - image_o[0]) = rho
    # cos(phi)*(x - PC[0]) + sin(phi)*(y - PC[1]) = tan(sigma)*PC[2]
    theta = SC.phi
    rho = (PC[1] - image_o[1])*np.cos(SC.phi) - (PC[0] - image_o[0])*np.sin(SC.phi) + np.tan(SC.sigma)*PC[2]
    while theta >= np.pi:
        theta -= np.pi
        rho *= -1.
    while theta < 0.0:
        theta += np.pi
        rho *= -1.
    return [rho, theta]

# l2data(Hough変換の rho方向2次微分)の極小値探索
def searchPeak(lx, # x-coordinates, 
#               ly, # y-coordinates,
               lyedata,   # Estimated errors of entries of ly
               l2data,    # 2nd derivatives of ly
               ibegin, iend,
               thred,     # "peak-height >= lyedata[i]*thred" is used as the criterion.
               peak_data): # Peak-searching result
    # H_FACTOR = -0.25*np.exp(0.5) # =-0.412180317675;
    # W_FACTOR = np.sqrt(2.0*np.log(2.0))
    H_FACTOR = -1.0 / (np.sqrt(3.0)*(np.exp(0.5)+2.0*np.exp(-1.0))) # = -0.24212836
    
    # Search for the first i with l2data[i] >= 0.0
    i = ibegin
    while i < iend:
        if l2data[i] >= 0.0:
            break
        i+=1
    
    peak_data.clear()
    while i < iend:
        # Search for the first i with l2data[i] < 0.0 
        if l2data[i] >= 0.0:
            i+=1
            continue
        
        i0 = i;	# l2data[i0-1] >= 0 and l2data[i0] < 0.
        Area_under0 = 0.0
        while i + 1 < iend:
            if l2data[i] > 0.0:
                break
            # Calculate the area surrounded by x-axe and the 2nd derivative under 0.
            Area_under0 += l2data[i] * (lx[i + 1] - lx[i - 1]) * 0.5
            i+=1
        
        # i > i0, l2data[i-1] <= 0 and l2data[i] > 0.
        if i + 1 >= iend:
            break
        
	    # Calculate x-values when 2nd derivative y=0.
        xp1 = (l2data[i0] * lx[i0 - 1] - l2data[i0 - 1] * lx[i0]) / (l2data[i0] - l2data[i0 - 1])
        xp2 = (l2data[i] * lx[i - 1] - l2data[i - 1] * lx[i]) / (l2data[i] - l2data[i - 1])
        dev2 = xp2 - xp1
        
        # Peak height obtained from the 2nd derivative.
        pkheight_area = H_FACTOR * dev2 * Area_under0
        
        j = i0
        while j < i:
            if pkheight_area < lyedata[j]*thred:
                break
            j += 1
        if j < i:
            i+=1
            continue
        
        # 2nd derivative is < 0 in (xp1, xp2) (almost equal to [l2data[i0], l2data[i])), peak-position, peak-height, fwhm
        pdata = PeakData()
        pdata.iinterval = [i0, i]
        pdata.dinterval = [xp1, xp2]
        pdata.peak_height = pkheight_area
        # pdata.fwhm = W_FACTOR * dev2      # Not used in this code.
        # print ([i0, i], [xp1, xp2], pkheight_area)
        peak_data.append(pdata)
        i+=1

# l1_H = image.shape[0] : 画像の高さ
# l2_W = image.shape[1] : 画像の幅 
def limit_rhos (l1_H, l2_W, rhos, theta):
    theta_rad = np.deg2rad (theta)
    # (rho, theta)は画像左上を原点とし、X軸を上, y軸を右向きに取ったとき以下の直線に対応:
    # -sin(theta)*(x - image_shape[0]//2) + cos(theta)*(y - image_shape[1]//2) = rho
    image_o = [l1_H//2, l2_W//2]
    if theta < 90.0:
        rho_min = -np.sin(theta_rad)*(l1_H - image_o[0]) + np.cos(theta_rad)*(     - image_o[1])
        rho_max = -np.sin(theta_rad)*(     - image_o[0]) + np.cos(theta_rad)*(l2_W - image_o[1])
    else:
        rho_min = -np.sin(theta_rad)*(l1_H - image_o[0]) + np.cos(theta_rad)*(l2_W - image_o[1])
        rho_max = -np.sin(theta_rad)*(     - image_o[0]) + np.cos(theta_rad)*(     - image_o[1])
    ibegin = bisect.bisect_right(rhos, rho_min)
    iend = bisect.bisect_left(rhos, rho_max)
    return ibegin, iend


# ArrayDeriv2(Hough変換の rho方向2次微分)を用いたバンド探索
def searchBand(rhos, # rho values
               thetas, # theta values
               ArrayDeriv2,    # 2nd derivatives of ly
               ArraySinogramErrors,   # Estimated errors of entries of ly
               image_shape, PC,
               BandKukans): # Output
    #import ebsd
    import params
    BandKukans.clear()
    # Returns an empty BandKukan if ArrayDeriv2 is empty.
    if len(ArrayDeriv2) < 1:
        return
    image_o = [image_shape[0]//2, image_shape[1]//2]
    rho_o = len(rhos)//2
    for t in range(len(ArrayDeriv2[0])):
        # X=[]
        ObsY0Error=[]
        Y2=[]
        for n in range(len(ArrayDeriv2)):
            # X.append(n) #x座標
            ObsY0Error.append(ArraySinogramErrors[n][t])
            Y2.append(ArrayDeriv2[n][t])
        Y2_ = np.array(Y2)*(-1)
        peak_data_min = [] # 2次微分の極小値 (=Hough変換の極大値, -)に関する情報を格納
        peak_data_max = [] # 2次微分の極大値 (=Hough変換の極小値, +)に関する情報を格納
        # ラドン変換の断面図におけるピークサーチ
        if params.Circle:
            searchPeak(rhos, ObsY0Error, Y2, 0, len(Y2), params.thred, peak_data_min)
            searchPeak(rhos, ObsY0Error, Y2_, 0, len(Y2_), params.thred, peak_data_max)
        else:
            ibegin, iend = limit_rhos (*image_shape, rhos, thetas[t])
            searchPeak(rhos, ObsY0Error, Y2, ibegin, iend, params.thred, peak_data_min)
            searchPeak(rhos, ObsY0Error, Y2_, ibegin, iend, params.thred, peak_data_max)
 
        # 2次微分の極大値のピーク(+)と極小値のピーク(-)の間にある零点の座標の検出
        # 各thetaにおける(+,-)型, (-,+)型 -> RhoOfBandEdges
        RhoOfBandEdges=[]
        j1=0
        j2=0
        for i in range(len(peak_data_min)):
            while j1 < len(peak_data_max) and peak_data_max[j1].iinterval[1] < peak_data_min[i].iinterval[0]:
                j1 = j1+1
            if j1 < len(peak_data_max) and peak_data_max[j1].iinterval[1] == peak_data_min[i].iinterval[0]:# +_-型
                edata = EdgeData()
                edata.left = peak_data_max[j1]
                edata.right = peak_data_min[i]
                edata.zero_rho = edata.left.dinterval[1]
                edata.isPlusMinus = True
                RhoOfBandEdges.append(edata)
            
            while j2 < len(peak_data_max) and peak_data_max[j2].iinterval[0] < peak_data_min[i].iinterval[1]:
                j2 = j2+1
            if j2 < len(peak_data_max) and peak_data_max[j2].iinterval[0] == peak_data_min[i].iinterval[1]:# -_+型
                edata = EdgeData()
                edata.left = peak_data_min[i]
                edata.right = peak_data_max[j2]
                edata.zero_rho = edata.right.dinterval[0]
                edata.isPlusMinus = False
                RhoOfBandEdges.append(edata)
        
        # 各thetaにおいて、零点の座標でソート
        RhoOfBandEdges.sort(key=lambda x: x.zero_rho)
    
        # バンドエッジの(rho, theta) -> rhotheta
        # 零点の間の区間に関する情報 -> BandKukans
        for m in range(len(RhoOfBandEdges)-1):
            edge1=RhoOfBandEdges[m]
            if not edge1.isPlusMinus: # +_- -_+型以外の排除
                continue
            bdata = None
            for m2 in range(m+1, len(RhoOfBandEdges)):
                edge2=RhoOfBandEdges[m2]
                x1=edge1.zero_rho 
                x2=edge2.zero_rho
                l=abs(x1-x2)
                if edge2.isPlusMinus: # +_- -_+型以外の排除
                    continue
                bdata2 = BandData()
                bdata2.setEdges(thetas[t], edge1.left.dinterval, edge2.right.dinterval)
                bdata2.setCenter(PC, image_o)
                if BAND_WIDTH_MIN > l:
                    continue
                if l > BAND_WIDTH_MAX:
                    break
                bm = mask.BandMask()
                bm.setParam(bdata2.putTheta(), bdata2.putEdgeRanges(), params.Circle)
                bdata2.convolution = bm.putConvolution(ArrayDeriv2, image_shape, rhos, rho_o, thetas)
                if params.MinCorrelation > bdata2.putConvolution():
                    continue
                if bdata is not None and bdata.putConvolution() > bdata2.putConvolution():
                    continue
                bdata = bdata2
            if bdata is not None:
                BandKukans.append(bdata)


def isEqualTheta(theta1, theta2, dtheta):
    dtheta2 = theta1 - theta2
    if abs( dtheta2 ) > dtheta and abs(dtheta2 + 180.0) > dtheta and abs(dtheta2 - 180.0) > dtheta:
         return False
    return True

# thetaの差 <= dthetaで、rhoの区間が重なるものの中で、peak_heightが最大のもののみをBandKukansに格納
def selectBands(dtheta, BandKukans): # Output
    copy = BandKukans.copy()
    BandKukans.clear()
    flags = [True]*len(copy)
    for i in range(len(copy)):
        if not flags[i]:
            continue
        band1 = copy[i]
        for j in range(len(copy)):
            if i == j:
                continue
            band2 = copy[j]
            if not isEqualTheta(band1.putTheta(), band2.putTheta(), dtheta):
                continue
            if abs(band1.putTheta() - band2.putTheta()) > 90.:
                if not isIntersected(band1.edge_rhos, [-band2.edge_rhos[1], -band2.edge_rhos[0]]):
                    continue
            else:
                if not isIntersected(band1.edge_rhos, band2.edge_rhos):
                    continue
            if band1.putConvolution() < band2.putConvolution():
                flags[i] = False
                break
            if band1.putConvolution() > band2.putConvolution():
                flags[j] = False
            elif( abs(band1.putTheta()   - band2.putTheta()  ) <= max(abs(band1.putTheta()),   abs(band2.putTheta())  )*1.0e-8 and 
                  abs(band1.edge_rhos[0] - band2.edge_rhos[0]) <= max(abs(band1.edge_rhos[0]), abs(band2.edge_rhos[0]))*1.0e-8 and 
                  abs(band1.edge_rhos[1] - band2.edge_rhos[1]) <= max(abs(band1.edge_rhos[1]), abs(band2.edge_rhos[1]))*1.0e-8 ):
                flags[j] = False
        if flags[i]:
            BandKukans.append(band1)


#与えられた2つの区間X,Yが共通部分を持つかの判定
def isIntersected(X, Y):
    if Y[1] < X[0] or X[1] < Y[0]:
        return False
    else:
        return True

# 昇順ソートされた配列arrの中でvalに最も近いインデックスを返す
def getNearestIndex(arr, val):
    i = bisect.bisect(arr, val)
    if i==0:
        return 0
    elif i==len(arr) or val - arr[i-1] <= arr[i] - val:
        return i-1
    else:
        return i

def getNearestZeroPoint(itheta, rho):
    global ArrayDeriv2
    min, max = 0, len(ArrayDeriv2)
    min_found, max_found = False, False
    for i in range(len(ArrayDeriv2) - 1):
        v1, v2 = ArrayDeriv2[i][itheta], ArrayDeriv2[i+1][itheta]
        if v1 == v2: continue
        pos = (v2*i - v1*(i+1)) / (v2 - v1)
        if v1>=0 and v2<=0 and i <= rho:
            min = pos
            min_found = True
        elif v1<=0 and v2>=0 and i > rho:
            max = pos
            max_found = True
            break
    if not min_found or not max_found: return None
    return [min - 0.5*len(ArrayDeriv2), max - 0.5*len(ArrayDeriv2)]

def getNearZeroPoint(itheta, rho, PlusMinus):
    global rhos, ArrayDeriv2
    edge=[0.,0.]
    i = getNearestIndex(rhos, rho)
    if PlusMinus:
        # +-型を探す
        if ArrayDeriv2[i][itheta] > 0.:
            while i<len(ArrayDeriv2) and ArrayDeriv2[i][itheta] > 0.: i+=1 
            if i>=len(ArrayDeriv2): return None
            i-=1
        else:
            while i>=0 and ArrayDeriv2[i][itheta] <= 0.: i-=1
            if i<0: return None
        # この時点で ArrayDeriv2[i][itheta] > 0 and ArrayDeriv2[i+1][itheta] <= 0.
        edge[1] = (ArrayDeriv2[i+1][itheta]*rhos[i] - ArrayDeriv2[i][itheta]*rhos[i+1]) / (ArrayDeriv2[i+1][itheta] - ArrayDeriv2[i][itheta])
        while i>=0 and ArrayDeriv2[i][itheta] > 0.: i-=1
        if i<0: return None
        # この時点で ArrayDeriv2[i] <= 0 and ArrayDeriv2[i+1] > 0.
        edge[0] = (ArrayDeriv2[i+1][itheta]*rhos[i] - ArrayDeriv2[i][itheta]*rhos[i+1]) / (ArrayDeriv2[i+1][itheta] - ArrayDeriv2[i][itheta])
    else:
        # -+型を探す
        if ArrayDeriv2[i][itheta] > 0.:
            while i>=0 and ArrayDeriv2[i][itheta] > 0.: i-=1 
            if i<0: return None
            i+=1
        else:
            while i<len(ArrayDeriv2) and ArrayDeriv2[i][itheta] <= 0.: i+=1
            if i>=len(ArrayDeriv2): return None
        # この時点で ArrayDeriv2[i-1][itheta] <= 0 and ArrayDeriv2[i][itheta] > 0.
        edge[0] = (ArrayDeriv2[i][itheta]*rhos[i-1] - ArrayDeriv2[i-1][itheta]*rhos[i]) / (ArrayDeriv2[i][itheta] - ArrayDeriv2[i-1][itheta])
        while i<len(ArrayDeriv2) and ArrayDeriv2[i][itheta] > 0.: i+=1
        if i>=len(ArrayDeriv2): return None
        # この時点で ArrayDeriv2[i-1][itheta] > 0. and ArrayDeriv2[i][itheta] <=0.
        edge[1] = (ArrayDeriv2[i][itheta]*rhos[i-1] - ArrayDeriv2[i-1][itheta]*rhos[i]) / (ArrayDeriv2[i][itheta] - ArrayDeriv2[i-1][itheta])
    return edge

#|
#| 出力処理用
#|

PC = None
thetas = None # degree
rhos = None
BandKukans = None
shape = None
ArraySinogramErrors = None
ArrayDeriv2 = None
name_pickle = 'out.pickle'
name_ArrayDeriv2 = 'out.2nd_derivative.tif'

def printAll():
    global PC
    global thetas
    global rhos
    global BandKukans
    global shape
    global ArraySinogramErrors
    with open(name_pickle, 'wb') as f:
        pickle.dump(PC, f)
        pickle.dump(thetas, f)
        pickle.dump(rhos, f)
        pickle.dump(BandKukans, f)
        pickle.dump(shape, f)
        pickle.dump(ArraySinogramErrors, f)
    printShapes(BandKukans, shape, 'out.shapes.json') # 可視化用テキストファイル出力
    printRhoTheta(BandKukans, 'out.rho_theta.txt') # θ, ρ_begin, ρ_end
    printSphericalCoordinates(BandKukans, 'data0.txt', flg = 0) # 抽出したバンドエッジ[rho_1, rho_2, theta]をテキスト出力
    printSphericalCoordinates(BandKukans, 'data1.txt', flg = 1)

def readResultsIfExists():
    global PC
    global thetas # degree
    global rhos
    global BandKukans
    global shape
    global ArraySinogramErrors
    global ArrayDeriv2
    try:
        if os.path.exists(name_pickle):
            with open(name_pickle, 'rb') as f:
                PC = pickle.load(f)
                thetas = pickle.load(f)
                rhos = pickle.load(f)
                BandKukans = pickle.load(f)
                shape = pickle.load(f)
                ArraySinogramErrors = pickle.load(f)
            ArrayDeriv2 = imread(name_ArrayDeriv2)
    except Exception as e:
        print(e, file=sys.stderr, flush=True)
        thetas = None
        rhos = None
        BandKukans = None
        shape = None
        ArraySinogramErrors = None
        ArrayDeriv2 = None
readResultsIfExists()


def printSphericalCoordinates(BandKukans, fname, flg = 0):
    with open(fname, 'w') as f:
        f.write('#	Use the band widths? (0: No, 1: Yes)   WaveLength(Angstrom)\n')
        f.write("{}	0.085885\n".format(flg))
        f.write('#	Phi(deg)		Sigma(deg)	Sigma_begin(deg)	Sigma_end(deg)\n')
        for band in BandKukans:
            f.write(f'{np.rad2deg(band.center_SC.phi):.5f}  {np.rad2deg(band.center_SC.sigma):.5f}  {np.rad2deg(band.center_SC.sigma-band.BraggAngle):.5f}  {np.rad2deg(band.center_SC.sigma+band.BraggAngle):.5f}\n')

# 配列ArrayRhoThetaに格納された各(rho, theta)=ArrayRhoTheta[i]に対応する、画像上の端点を算出．
# 座標(rho, theta)は実数値とする。
# 上記のgetLineForDisplay関数を使用．
def getLinesForDisplay(image_shape, ArrayRhoTheta, ArrayLines):
    ArrayLines.clear()
    for rt in ArrayRhoTheta:
        LineX = []
        LineY = []
        getLineForDisplay(image_shape, rt[0], np.deg2rad(rt[1]), LineX, LineY)
        # 端点はちょうど２つのはず
        # ちょうど画像の角を通る）
        if len(LineX) != 2: continue
        ArrayLines.append([LineX, LineY])

def printShapes(BandKukans, shape, fname):
    BandEdges_rhotheta  = []
    BandCenters_rhotheta = []
    for band in BandKukans:
        BandEdges_rhotheta.append([band.edge_rhos[0], band.center_rt[1]])
        BandEdges_rhotheta.append([band.edge_rhos[1], band.center_rt[1]])
        BandCenters_rhotheta.append(band.center_rt)

    edges = []
    centers = []
    getLinesForDisplay(shape, BandEdges_rhotheta, edges)
    getLinesForDisplay(shape, BandCenters_rhotheta, centers)

    with open(fname, 'w') as f:
        f.write('{\n')
        f.write('  "band-edge": [\n')
        for Line in edges:
            f.write(f'    {{"from": [{Line[0][0]},{Line[1][0]}], "to": [{Line[0][1]},{Line[1][1]}]}}{"" if Line==edges[-1] else ","}\n')
        f.write('  ],\n')
        f.write('  "band-center": [\n')
        for Line in centers:
            f.write(f'    {{"from": [{Line[0][0]},{Line[1][0]}], "to": [{Line[0][1]},{Line[1][1]}]}}{"" if Line==centers[-1] else ","}\n')
        f.write('  ],\n')
        f.write('  "band-edge-2nd": [\n')
        for i in range(0, len(BandEdges_rhotheta), 2):
            f.write(f'    {{"from": [{BandEdges_rhotheta[i][1]},{BandEdges_rhotheta[i][0]}], "to": [{BandEdges_rhotheta[i+1][1]},{BandEdges_rhotheta[i+1][0]}]}}{"" if i==len(BandEdges_rhotheta)-2 else ","}\n')
        f.write('  ]\n')
        f.write('}\n')

def printRhoTheta(BandKukans, fname):
    if BandKukans is None: return
    with open(fname, 'w') as f:
        f.write('#  theta  rho_center  rho_begin  rho_end  score\n')
        for band in BandKukans:
            f.write(f'{band.center_rt[1]} {band.center_rt[0]} {band.edge_rhos[0]} {band.edge_rhos[1]} {band.convolution}\n')

# 計算結果を画面上で確認するための関数群
#[r,t] in rhotheta に対応する直線と元画像を重ねた画像を表示
"""def PlotLines(image, rhotheta, title):
    ArrayLines = []
    getLinesForDisplay(image.shape, rhotheta, ArrayLines)
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(1, 1, 1)
    plt.title('EBDS image and ' + title)
    for Line in ArrayLines:
        plt.plot(Line[0], Line[1], c='y')
    ax.imshow(image, cmap=plt.cm.Greys_r)

# バンドに対応する線分と二次微分を重ねて表示
def PlotLines2(ArrayDeriv2, rhotheta, title):
    fig = plt.figure(figsize=(4.5,4.5))
    ax = fig.add_subplot(1,1,1)
    plt.title('Second derivative and ' + title)
    for i in range(0, len(rhotheta), 2):
        plt.plot([rhotheta[i][1], rhotheta[i+1][1]], [rhotheta[i][0], rhotheta[i+1][0]], c='y')
    ax.imshow(ArrayDeriv2.tolist(), cmap=plt.cm.Greys_r,
           extent=(0, 180, len(ArrayDeriv2)*0.5, -len(ArrayDeriv2)*0.5), aspect='auto')

def plotLinesAndOriginal(image, rhotheta):
    ArrayLines = []
    getLinesForDisplay(image.shape, rhotheta, ArrayLines)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.5))
    for Line in ArrayLines:
        ax1.plot(Line[0], Line[1], c='y')
    ax1.imshow(image, cmap=plt.cm.Greys_r)
    ax2.imshow(image, cmap=plt.cm.Greys_r)"""

def setZeroOutsideCircle (image):
    shape_min = min(image.shape)
    radius = shape_min // 2
    img_shape = np.array(image.shape)
    coords = np.array(np.ogrid[:image.shape[0], :image.shape[1]],
                        dtype=object)
    dist = ((coords - img_shape // 2) ** 2).sum(0)
    outside_reconstruction_circle = dist > radius ** 2
    if np.any (image[outside_reconstruction_circle]):
        image[outside_reconstruction_circle] = 0.0
    return image, ~outside_reconstruction_circle

def calcSigma (image, mask):
    mean = image[mask].mean()
    sigma = np.mean (np.abs (image[mask] - mean))
    return sigma

#|
#| バンド抽出計算
#|

def run():
    global PC      # project centerの座標（3次元ベクトル, スケール変換後）
    global Circle  # EBSD画像が円かどうか
    global rhos # Hough変換のrho座標を格納した配列
    global thetas # Hough変換のtheta座標を格納した配列
    global BandKukans
    global shape #　EBSD画像の縦横サイズ(スケール変換後)
    global ArrayDeriv2 # Hough変換の rho方向2次微分値を格納した配列
    global ArraySinogramErrors  # Hough変換の誤差見積もり値を格納した配列
    global BAND_WIDTH_MIN
    global BAND_WIDTH_MAX

    try:
        import file
        #import ebsd
        import params
        importlib.reload (file)   # file.pyの読み込み
        importlib.reload (params) # params.pyの読み込み
        
        # 入力ファイル指定
        filename = file.path     # EBSD画像ファイルの　path 
        PC0 = params.PC0         # 下式でproject centerの座標（3次元ベクトル, スケール変換前）は求められるとする:
        Circle = params.Circle   # True: EBSD画像は円, False: 四角
        #print (filename)
        
        time_st = time.time()
        # 画像のリスケール
        print('Rescale image...', flush=True)
        time_st_pr = time.time()
        image = imread(filename, as_gray=True) # 画像の読込み
        RescaleParam = params.RescaleParam / max(image.shape) # 画像のスケールを縮小するパラメータ
        image = rescale(image, scale=RescaleParam, mode='reflect') # 画像のスケールを変更
        mask_circle = np.ones_like (image).astype (np.bool_)
        if Circle:
            image, mask_circle = setZeroOutsideCircle (image)
        # sigma(EBSD画像の標準偏差)を推定
        print('Calculate error...', flush=True)
        sigma = calcSigma (image, mask_circle)
        # L=[]
        # for n in range(len(image)):
        #     L.extend(image[n])
        # m = sum(L)/len(L) #画像の平均
        # print(f'  Mean of input image = {m}', flush=True)
        # LL = abs(L - m)
        # sigma = sum(LL)/len(LL) #画像の標準偏差 => 絶対値の平均（外れ値の影響軽減のため）
        print(f'  Mean error of input image = {sigma}', flush=True)
        print(f'  Error estimated by using the input threshold = {sigma*params.thred}', flush=True)

        # band width min maxを画像サイズから再計算
        BAND_WIDTH_MIN = params.BAND_WIDTH_MIN * max(image.shape)
        BAND_WIDTH_MAX = params.BAND_WIDTH_MAX * max(image.shape)
        print ('updated band width min : {:.2f} px, max : {:.2f} px'.format(BAND_WIDTH_MIN, BAND_WIDTH_MAX))

        shape = image.shape
        imsave('out.rescaled.png', exposure.rescale_intensity(image, in_range='image', out_range='uint8').astype(np.uint8))
        PC = [PC0[1]*image.shape[0], PC0[0]*image.shape[1], PC0[2]*image.shape[0]] # PCの座標（スケール変換後）
        print(f'  Size: {shape} px', flush=True)
        print(f'  Projection center: {PC} px', flush=True)
        print(f'  EBSD image is a circle?: {Circle}', flush=True)
        print(f'  done', flush=True)
        
        # ラドン変換
        print('Radon transform...', flush=True)
        thetas = np.linspace(0., 180., max(image.shape), endpoint=False) # thetaのbin幅 = 180./画像の横幅・縦幅の大きい方, とする。
        sinogram = radon(image, theta=thetas, circle=Circle)              # ラドン変換
        imsave('out.radon.png', exposure.rescale_intensity(sinogram, in_range='image', out_range='uint8').astype(np.uint8))
        print(f'  Size: {sinogram.shape}', flush=True)
        print('  done', flush=True)
        
        # ラドン変換の誤差値(各(rho, sigma)に対応する(画像内線分の長さ)^{1/2}*sigma)を求める
        #　線分の長さはスケール変換後の値。
        ArraySinogramErrors = []
        calRadonGosa(thetas, image, sinogram, sigma    , Circle, ArraySinogramErrors    )
        print('  done', flush=True)
        
        # 微分値の計算
        print('Calculate derivatives...', flush=True)
        ArraySmth = np.zeros((sinogram.shape[0], sinogram.shape[1])) # 2D arrays for storing (rho, theta)
        ArrayDeriv1 = np.zeros((sinogram.shape[0], sinogram.shape[1]))
        ArrayDeriv2 = np.zeros((sinogram.shape[0], sinogram.shape[1]))
        rho_o = sinogram.shape[0]//2
        rhos = np.arange(-rho_o, sinogram.shape[0]-rho_o) # sinogram.shape[0] = rhoのbin数.
        for k in np.arange(sinogram.shape[1]): # sinogram.shape[1] = thetaのbin数.
            Y = [] # sinogramの列ベクトル(theta=const.)を格納
            for n in np.arange(sinogram.shape[0]):
                Y.append(sinogram[n][k])
            a=0
            b=sinogram.shape[0]-1
            while a < len(sinogram) and sinogram[a][k] <= 0.:
                a+=1
            while b >= a and sinogram[b][k] <= 0.:
                b-=1
            for n in np.arange(max(a, params.num_points), 
                               min(b+1, sinogram.shape[0]-params.num_points)):
                coef3 = putCoef(rhos,Y,n,params.num_points,params.deg)
                ArrayDeriv2[n][k] = coef3[1]*2.0
                ArrayDeriv1[n][k] = coef3[2]
                ArraySmth[n][k] = coef3[3]
        time_ed_pr = time.time()
        #print ('time for preprocess {} sec'.format (time_ed_pr - time_st_pr))
        imsave('out.1st_derivative.png', exposure.rescale_intensity(ArrayDeriv1, in_range='image', out_range='uint8').astype(np.uint8))
        imsave('out.2nd_derivative.png', exposure.rescale_intensity(ArrayDeriv2, in_range='image', out_range='uint8').astype(np.uint8))
        imsave(name_ArrayDeriv2, ArrayDeriv2)
        print('  done', flush=True)
        
		# 計算結果を画面上で確認するため
        # 1次微分, 2次微分の画像作成
        """fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.5))
        ax1.set_title('First derivative')
        ax1.set_xlabel("Projection angle (deg)")
        ax1.set_ylabel("Projection position (pixels)")
        ax1.imshow(ArrayDeriv1.tolist(), cmap=plt.cm.Greys_r,
            extent=(0, 180, len(ArrayDeriv1)*0.5, -len(ArrayDeriv1)*0.5), aspect='auto')
        ax2.set_title('Second derivative')
        ax2.set_xlabel("Projection angle (deg)")
        ax2.set_ylabel("Projection position (pixels)")
        ax2.imshow(ArrayDeriv2.tolist(), cmap=plt.cm.Greys_r,
            extent=(0, 180, len(ArrayDeriv2)*0.5, -len(ArrayDeriv2)*0.5), aspect='auto')
        fig.tight_layout()"""
        #plt.show()
        
        # バンド抽出
        print('Search bands...', flush=True)
        BandKukans = []
        time_st_bs = time.time()
        searchBand(rhos, thetas, ArrayDeriv2, ArraySinogramErrors    , 
                   image.shape, PC, BandKukans    )
        print(f'  bands: {len(BandKukans)} => ', end='', flush=True)
        # thetaの差がdtheta(degree)以下で、rhoの区間が重なるバンドのうち　convolutionの値が最大値となるバンドのリスト
        selectBands(params.dtheta, BandKukans)
        
        print(len(BandKukans), flush=True)
        BandKukans.sort(key=lambda b: b.putConvolution(), reverse=True)
        time_ed_bs = time.time()
        #print ('\ntime for band search {} sec'.format(time_ed_bs - time_st_bs))
        # iter = 0
        # while iter < 1 and params.GenerateBandsFromBands and len(BandKukans) < params.NumberOfBands:
        #    iter += 1
        #    addBandsFrom4Bands()
        
        # BandDataの全てのパラメータをセット
        BandEdges_rhotheta  = [] # PlotLinesのための配列
        BandCenters_rhotheta = [] # PlotLinesのための配列
        for band in BandKukans:
            BandEdges_rhotheta.append([band.edge_rhos[0], band.center_rt[1]])
            BandEdges_rhotheta.append([band.edge_rhos[1], band.center_rt[1]])
            BandCenters_rhotheta.append(band.center_rt)
        #PlotLines2(ArrayDeriv2, BandEdges_rhotheta, 'Band edges') # Hough変換の2次微分とバンド区間に対応する線を重ねたものを表示
        # PlotLines(image, BandCenters_rhotheta, 'Band centers')  # EBSD画像とバンドセンターに対応する直線を重ねたものを表示
        #plotLinesAndOriginal(image, BandEdges_rhotheta)           # EBSD画像とバンドエッジに対応する直線を重ねたものを表示
        #plotLinesAndOriginal(image, BandCenters_rhotheta)
        # 出力した画像をpdfで保存
        #pdf = PdfPages(filename + '.out.pdf')
        #fignums = plt.get_fignums()
        #for fignum in fignums:
        #    plt.figure(fignum)
        #    pdf.savefig()
        #pdf.close()
        
        # 結果を出力
        printAll()
        time_ed = time.time ()
        #print ('all process time {}'.format(time_ed - time_st))

    except Exception as e:
        print(e, file=sys.stderr, flush=True)
    
    finally:
        print("--- done ---", flush=True)

#|
#| ex) removeBands([2,0,4])  # 0,2,4番目のバンドを消去する
#|
def removeBands(indices):
    global BandKukans
    if BandKukans is None: return
    for i in sorted(indices, reverse=True):
        del BandKukans[i]
    printAll()

#|
#| BandKukans[i]のバンド中心rhoを変更する(バンドセンター以外の情報は変更なし)
#|
def editBandCenter(rho, i):
    global PC, BandKukans, shape, ArrayDeriv2
    if ArrayDeriv2 is None: return
    BandKukans[i].setCenter(PC, shape, rho)
    printAll()

#|
#| BandKukans内にバンド(バンドセンター: theta, rho)が存在すればそれを返す。なければNone。
#|
def find(theta, rho, BandKukans):
    import params
    for b in BandKukans:
        match_theta = isEqualTheta(theta, b.putTheta(), params.dtheta)
        if abs(theta - b.putTheta()) > 90.:
            match_rho = b.edge_rhos[0]-0.001 <= -rho and -rho <= b.edge_rhos[1]+0.001
        else:
            match_rho = b.edge_rhos[0]-0.001 <= rho and rho <= b.edge_rhos[1]+0.001
        # 0.001だけ広げる理由：グラフ上で画像ではなくバンド線をクリックするとバンドの境界値が取得されてしまう
        if match_rho and match_theta:
            return b
    return None

#|
#| BandKukans内に同じバンドが存在すればそれを返す。なければNone。
#|
def findBand(band, BandKukans):
    import params
    for b in BandKukans:
        match_theta = isEqualTheta(band.putTheta(), b.putTheta(), params.dtheta)
        ranges1 = band.putEdgeRanges()
        ranges2 = b.putEdgeRanges()
        if abs(band.putTheta() - b.putTheta()) > 90.:
            match_rho = (isIntersected(ranges1[0], [-ranges2[1][1], -ranges2[1][0]]) and
                         isIntersected(ranges1[1], [-ranges2[0][1], -ranges2[0][0]]))
        else:
            match_rho = isIntersected(ranges1[0], ranges2[0]) and isIntersected(ranges1[1], ranges2[1])
        if match_rho and match_theta:
            return b
    return None

def findAll(band, BandKukans, dtheta):
    ans = []
    for i in range(len(BandKukans)):
        match_theta =  isEqualTheta(band.putTheta(), BandKukans[i].putTheta(), dtheta)
        ranges1 = band.putEdgeRanges()
        ranges2 = BandKukans[i].putEdgeRanges()
        if abs(band.putTheta() - BandKukans[i].putTheta()) > 90.:
            match_rho = (isIntersected(ranges1[0], [-ranges2[1][1], -ranges2[1][0]]) and 
                         isIntersected(ranges1[1], [-ranges2[0][1], -ranges2[0][0]]))
        else:
            match_rho = isIntersected(ranges1[0], ranges2[0]) and isIntersected(ranges1[1], ranges2[1])
        if match_rho and match_theta:
            ans.append(i)
    return ans

# バンド中心（角度θ、中心ρ）からバンドを生成する
def getBand_theta_rho(theta0, rho0, BAND_WIDTH_MIN, BAND_WIDTH_MAX):
    global PC, Circle, rhos, thetas, BandKukans, shape, ArrayDeriv2
    if ArrayDeriv2 is None: return
    rho_o = len(rhos)//2
    itheta = getNearestIndex(thetas, theta0)
    irho_begin1 = bisect.bisect(rhos, rho0 - BAND_WIDTH_MAX*0.667)
    irho_end1   = bisect.bisect(rhos, rho0 - BAND_WIDTH_MIN*0.333)
    irho_begin2 = bisect.bisect(rhos, rho0 + BAND_WIDTH_MIN*0.333)
    irho_end2   = bisect.bisect(rhos, rho0 + BAND_WIDTH_MAX*0.667)
    ObsY0Error = []
    Y2=[]
    for n in range(len(ArrayDeriv2)):
        Y2.append(-ArrayDeriv2[n][itheta])
        ObsY0Error.append(ArraySinogramErrors[n][itheta])
    peak_data_left = [] # 2次微分の極大値に関する情報を格納
    peak_data_right = [] # 2次微分の極大値に関する情報を格納
    searchPeak(rhos, ObsY0Error, Y2, irho_begin1, irho_end1, 0.001, peak_data_left)
    searchPeak(rhos, ObsY0Error, Y2, irho_begin2, irho_end2, 0.001, peak_data_right)
    bdata = None
    for peak_left in peak_data_left:
        for peak_right in peak_data_right:
            drho = peak_right.dinterval[0] - peak_left.dinterval[1]
            if drho < BAND_WIDTH_MIN or BAND_WIDTH_MAX < drho:
                continue
            if rho0 < peak_left.dinterval[1] + drho*0.35 or peak_right.dinterval[0] - drho*0.35 < rho0:
                continue
            bdata2 = BandData()
            bdata2.setEdges(thetas[itheta], peak_left.dinterval, peak_right.dinterval)
            bdata2.setCenter(PC, shape)
            bm = mask.BandMask()
            bm.setParam(bdata2.putTheta(), bdata2.putEdgeRanges(), Circle)
            bdata2.convolution = bm.putConvolution(ArrayDeriv2, shape, rhos, rho_o, thetas)
            if bdata is not None and bdata.putConvolution() > bdata2.putConvolution():
                continue
            bdata = bdata2
    return bdata

#|
#| バンドエッジ（theta, rho1, rho2）からバンドを生成する
#|
def getBand_theta_rhos(theta, rho1, rho2):
    global PC, Circle, rhos, thetas, BandKukans, shape, ArrayDeriv2
    if ArrayDeriv2 is None: return
    # BandKukansの中からバンドの候補を探す
    itheta = getNearestIndex(thetas, theta)
    edge_range1 = getNearZeroPoint(itheta, rho1, True)
    if edge_range1 is None:
        return None
    edge_range2 = getNearZeroPoint(itheta, rho2, False)
    if edge_range2 is None:
        return None
    band = BandData()
    band.setEdges(thetas[itheta], edge_range1, edge_range2)
    band.setCenter(PC, shape)
    bm = mask.BandMask()
    bm.setParam(band.putTheta(), band.putEdgeRanges(), Circle)
    rho_o = len(rhos)//2
    band.convolution = bm.putConvolution(ArrayDeriv2, shape, rhos, rho_o, thetas)
    return band

#|
#| バンド（角度θ、中心ρ）をBandKukansに追加する。
#|
def addBand_theta_rho(theta, rho):
    global BandKukans
    # 同じバンドが既に存在する場合は何もしない
    if find(theta, rho, BandKukans) is not None:
        print('Failed: band already exists', file=sys.stderr, flush=True)
        return
    # バンドを取得して追加する
    band = getBand_theta_rho(theta, rho, BAND_WIDTH_MIN, BAND_WIDTH_MAX)
    if band is None:
        print('Failed: no band found', file=sys.stderr, flush=True)
        return
    BandKukans.append(band)
    BandKukans.sort(key=lambda b: b.putConvolution(), reverse=True)
    print(f'ADDED: (θ, ρ_cener, ρ_begin, ρ_end) = ({band.center_rt[1]:.4f} {band.center_rt[0]:.4f} {band.edge_rhos[0]:.4f} {band.edge_rhos[1]:.4f})', flush=True)
    printAll()

#|
#| バンド（角度targetTheta(degree)、境界[rhomin,rhomax]）をBandKukansに追加する。
#|
def addBand_theta_edges(targetTheta, rhomin, rhomax):
    global BandKukans
    # バンドを取得
    band = getBand_theta_rhos(targetTheta, rhomin, rhomax)
    if band is None:
        print('Failed: no band found', file=sys.stderr, flush=True)
        return
    # 同じバンドが既に存在する場合は何もしない
    if findBand(band, BandKukans) is not None:
        print('Failed: band already exists', file=sys.stderr, flush=True)
        return
    # バンドを追加する
    BandKukans.append(band)
    BandKukans.sort(key=lambda b: b.putConvolution(), reverse=True)
    print(f'ADDED: (θ, ρ_cener, ρ_begin, ρ_end) = ({band.center_rt[1]:.4f} {band.center_rt[0]:.4f} {band.edge_rhos[0]:.4f} {band.edge_rhos[1]:.4f})', flush=True)
    printAll()

#|
#| 2つのバンドセンターの交点を返す
#|
def getCrossing(band1, band2):
    global Circle, shape
    import params
    rho1, theta1 = band1.center_rt
    rho2, theta2 = band2.center_rt
    # 角度が同じだと交点がないのでNoneを返す
    if isEqualTheta(theta1, theta2, params.dtheta):
        return None
    # 交点の計算
    Line1X, Line1Y = [], []
    Line2X, Line2Y = [], []
    mask.getLine(shape, rho1, np.deg2rad(theta1), Line1X, Line1Y, Circle) # 端点の座標をp0, p1とする
    mask.getLine(shape, rho2, np.deg2rad(theta2), Line2X, Line2Y, Circle) # 端点の座標をq0, q1とする
    # M = (-p1+p0, q1-q0), b=p0-q0
    # x=M^{-1}*bとすると、p0 + x[0]*(p1-p0) = q0 + x[1]*(q1-q0)
    M = [[-Line1X[1]+Line1X[0], Line2X[1]-Line2X[0]], [-Line1Y[1]+Line1Y[0], Line2Y[1]-Line2Y[0]]]
    det = M[0][0]*M[1][1] - M[0][1]*M[1][0]
    c = (M[1][1] * (Line1X[0]-Line2X[0]) - M[0][1] * (Line1Y[0]-Line2Y[0])) / det # c = x[0]
    # pos = p0 + c(p1-p0)
    pos = [Line1X[0] + c*(Line1X[1]-Line1X[0]), Line1Y[0] + c*(Line1Y[1]-Line1Y[0])]
    return pos


#|
#| 4つのバンドの交点を通るバンドを追加する
#|
def addBandsFrom4BandsIn(BandKukans, BAND_WIDTH_MIN, BAND_WIDTH_MAX, MinCorrelation, newBands):
    global shape
    size = len(BandKukans)
    image_o = [shape[0]//2, shape[1]//2]
    newBands.clear()
    for i in range(0, size):
        for j in range(i+1, size):
            for k in range(j+1, size):
                for l in range(k+1, size):
                    # バンドの交点の計算
                    crossings = [
                        getCrossing(BandKukans[i], BandKukans[j]),
                        getCrossing(BandKukans[i], BandKukans[k]),
                        getCrossing(BandKukans[i], BandKukans[l]),
                        getCrossing(BandKukans[j], BandKukans[k]),
                        getCrossing(BandKukans[j], BandKukans[l]),
                        getCrossing(BandKukans[k], BandKukans[l]),
                    ]
                    # 交点を通る直線を取得
                    rts = [
                        mask.getRhoThetaFromCrossings(image_o, crossings[0], crossings[5]),
                        mask.getRhoThetaFromCrossings(image_o, crossings[1], crossings[4]),
                        mask.getRhoThetaFromCrossings(image_o, crossings[2], crossings[3]),
                    ]
                    # バンドの追加
                    for rho, theta in rts:
                        if theta is None: continue
                        # 既に存在している場合は何もしない
                        if find(np.rad2deg(theta), rho, BandKukans) is not None:
                            # print("IGNORED (already exists)", flush=True)
                            continue
                        # バンドが画像内を通らない場合は何もしない
                        # LineX, LineY = [], []
                        # mask.getLine(shape, rho, theta, LineX, LineY, False)
                        # if len(LineX) != 2:
                            # print("IGNORED (outside the image)", flush=True)
                            # continue
                        # バンドを追加
                        band = getBand_theta_rho(np.rad2deg(theta), rho, BAND_WIDTH_MIN, BAND_WIDTH_MAX)
                        if band is None or band.putConvolution() < MinCorrelation:
                            # print("IGNORED (not a band)", flush=True)
                            continue
                        else:
                            newBands.append(band)
                            print(f"  (θ, ρ) = ({np.rad2deg(theta):.3f}, {rho:.3f}) => ", end='', flush=True)
                            print("ADDED", flush=True)


def addBandsFrom4Bands():
    global BandKukans
    import params
    newBands = []
    addBandsFrom4BandsIn(BandKukans, BAND_WIDTH_MIN, BAND_WIDTH_MAX, params.MinCorrelation, newBands)
    BandKukans.extend(newBands)
    print(f'  bands: {len(BandKukans)} => ', end='', flush=True)
    selectBands(params.dtheta, BandKukans)
    print(len(BandKukans), flush=True)
    BandKukans.sort(key=lambda b: b.putConvolution(), reverse=True)
    printAll()
    print("--- done ---", flush=True)
    
if __name__ == '__main__':
    run()
