import numpy as np
import bisect
import importlib
#from matplotlib import pyplot as plt

# (rho, theta)に対応する直線の、画像上の端点を算出．
# (rho, theta)は画像を配列とみたときの座標系、すなわち、
# 画像左上を原点、X軸を下, y軸を右向きに取ったとき、以下の直線に対応:
# -sin(theta)*(x - image_shape[0]//2) + cos(theta)*(y - image_shape[1]//2) = rho
# すなわち indexが(image_shape[0]//2, image_shape[1]//2)に等しいpixelが原点に対応する.
def getLine (image_shape, rho, theta, lineX, lineY, circle : bool):
    if circle:
        getLineCircle (image_shape, rho, theta, lineX, lineY)
    else:
        getLineBox (image_shape, rho, theta, lineX, lineY)

def getLineBox(image_shape, rho, theta, LineX, LineY):
    LineX.clear()
    LineY.clear()
    image_o = [image_shape[0]//2,image_shape[1]//2]
    # 画像のエッジ y = 0, y = image.shape[1]との交点を求める.
    if theta != 0.0:
        y = np.array([0.0, image_shape[1]])
        x = image_o[0] + (y - image_o[1])/np.tan(theta) - rho/np.sin(theta) 
        for i in np.arange(0,2):
            if 0 <= x[i] and x[i] <= image_shape[0]:
                LineX.append(x[i])
                LineY.append(y[i])
    
    # 画像のエッジ x = 0, x = image.shape[0]との交点を求める
    if theta != np.pi*0.5:
        x = np.array([0.0, image_shape[0]])
        y = image_o[1] + (x - image_o[0])*np.tan(theta) + rho/np.cos(theta) 
        for i in np.arange(0,2):
            if 0 <= y[i] and y[i] <= image_shape[1]:
                LineX.append(x[i])
                LineY.append(y[i])

# rho : EBSD画像の中心-直線の距離
# theta : 法線ベクトルの角度 (rad)
def getLineCircle (image_shape, rho, theta, LineX, LineY):
    radius = min (image_shape) // 2 # skimage.radonでの半径計算式
    LineX.clear ()
    LineY.clear ()
    
    # 直線の法線と、直線と円の交点の角度αを求める。
    # rhoが微小値だけradiusより大きくなる場合がある。
    #alpha = 0.
    #if abs (rho) < radius:
    #    alpha = np.arccos (abs (rho) / radius)
    rho = np.sign (rho) * min (abs (rho), radius)
    alpha = np.arccos (rho / radius)
    
    # 2つの交点（直線と円周）の座標の計算
    # 画像中心座標を(xo,yo)とすると下記のようになる。
    # (x1, y1) = (- radius * sin(theta+α) + xo, radius * cos(theta+α) + yo)
    # (x2, y2) = (- radius * sin(theta-α) + xo, radius * cos(theta-α) + yo)
    # In this case,
    # -sin(theta)*(x - xo) + cos(theta)*(y - yo) = radius * cos(α) = rho
    LineX += [-radius * np.sin (theta + alpha) + image_shape[0] // 2,
              -radius * np.sin (theta - alpha) + image_shape[0] // 2]
    LineY += [ radius * np.cos (theta + alpha) + image_shape[1] // 2,
               radius * np.cos (theta - alpha) + image_shape[1] // 2]

#|
#| p1, p2を通る直線に対応する(rho, theta)を返す。
#|
def getRhoThetaFromCrossings(image_o, p1, p2):
    if p1 is None or p2 is None: 
        return [None, None]
    d = [p1[0]-p2[0], p1[1]-p2[1]]
    # (p1, p2)を通る直線: (p2[1] - p1[1])*(x - p1[0]) - (p2[0] - p1[0])*(y - p1[1]) = 0.
    # (rho, theta)を通る直線 -sin(theta)*(x - image_o[0]) + cos(theta)*(y - image_o[1]) = rho.
    # arctan2(y,x) := arctan(y/x). したがって、
    if d[1] == 0.0:
        theta = 0.0
    else:
        theta = np.arctan2(d[1], d[0]) if d[1] > 0.0 else np.arctan2(-d[1],-d[0])
    rho = -np.sin(theta) * (p1[0] - image_o[0]) + np.cos(theta) * (p1[1] - image_o[1])
    return [rho, theta]


class BandMask:
    def __init__(self):
        self.theta0 = 0.0 # radian
        self.sigma0 = 0.0 # px
        self.rho0 = [0.0, 0.0] # rho0[0] < rho1[0], px
        self.rho0_wide = [0.0, 0.0] # rho0_wide[0] < rho0[0] < rho1[0] < rho1_wide[0], px
        self.Circle    = False

    def setParam(self, theta0, rho0_ranges, circle):
        self.theta0 = np.deg2rad(theta0)
        drho0 = ((rho0_ranges[0][1] - rho0_ranges[0][0])+(rho0_ranges[1][1] - rho0_ranges[1][0]))*0.5
        self.rho0 = [rho0_ranges[0][1], rho0_ranges[1][0]]
        self.rho0_wide = [self.rho0[0]-drho0, self.rho0[1]+drho0]
        if (self.rho0_wide[0]*2.+self.rho0_wide[1])/3. < self.rho0[0]:
            self.rho0[0] = (self.rho0_wide[0]*2.+self.rho0_wide[1])/3. 
        if (self.rho0_wide[0]+self.rho0_wide[1]*2.)/3. > self.rho0[1]:
            self.rho0[1] = (self.rho0_wide[0]+self.rho0_wide[1]*2.)/3. 
        self.sigma0 = drho0/np.sqrt(3.)
        self.Circle = circle

    def MaskFunc1(self, rho, rho0_arg, sigma0_arg):
        diff = (rho - rho0_arg)/sigma0_arg
        return -np.exp(-0.5*(diff**2))*diff*(diff**2 - 3.)/(sigma0_arg**2)

    def MaskFunc2(self, rho, rho0_arg, sigma0_arg, theta, z):
        dtheta = theta - self.theta0
        cosd = np.cos(dtheta)
        sind = np.sin(dtheta)
        diff = (rho*cosd + z*sind - rho0_arg)/sigma0_arg
        return np.exp(-0.5*(diff**2))*(diff**2 - 1.)*(cosd**2)/sigma0_arg

    def MaskFunc3(self, rho, rho0_arg, sigma0_arg, theta, MinZ, MaxZ):
        if abs(theta - self.theta0) < 1.0e-8:
            return self.MaskFunc1(rho, rho0_arg, sigma0_arg)*(MaxZ-MinZ) # MaxZ - MinZ = 線分の長さ
        else:
            dtheta = theta - self.theta0
            sind = np.sin(dtheta)
            return (self.MaskFunc2(rho, rho0_arg, sigma0_arg, theta, MaxZ) - self.MaskFunc2(rho, rho0_arg, sigma0_arg, theta, MinZ))/sind

    def calMaskValue(self, image_shape, rho, theta):
        LineX=[]
        LineY=[]
        # 線分 -sin(theta)*(x - image_shape[0]//2) + cos(theta)*(y - image_shape[1]//2) = rhoの端点
        getLine(image_shape, rho, theta, LineX, LineY, self.Circle)
        # 線分上の座標
        x0 = image_shape[0]//2 - rho*np.sin(theta)
        y0 = image_shape[1]//2 + rho*np.cos(theta)
        # 線分は以下でパラメトライズされる。線分と画像と交わるzの範囲を求める。　
        # (x0 + z*cos(theta), y0 + z*sin(theta))
        # z=0の座標は画像上にあることを仮定。
        MinZ = 0.
        MaxZ = 0.
        for i in range(len(LineX)):
            z = (LineX[i]-x0)*np.cos(theta) + (LineY[i]-y0)*np.sin(theta) # zは端点(LineX[i], LineY[i])に対応
            MaxZ = max(MaxZ, z)
            MinZ = min(MinZ, z)
        return self.MaskFunc3(rho, self.rho0[1], self.sigma0, theta, MinZ, MaxZ) - self.MaskFunc3(rho, self.rho0[0], self.sigma0, theta, MinZ, MaxZ)

    # thetas and theta are in degrees
    def calRescaledValue(self, arr, rs, r_o, thetas, r, theta):
        if theta < 0.0:
            return self.calRescaledValue(arr, rs, r_o, thetas, -r, theta+180.0)
        elif theta >= 180.0:
            return self.calRescaledValue(arr, rs, r_o, thetas, -r, theta-180.0)
        
        ir1 = bisect.bisect(rs, r)
        if ir1 == 0:
            ir0 = ir1 
            wr = [1., 1.]
        elif ir1 == len(rs):
            ir0 = ir1-1 
            ir1 = ir0 
            wr = [1., 1.]
        else:
            ir0 = ir1-1 
            wr = [rs[ir1]-r, r-rs[ir0]]
        
        it1 = bisect.bisect(thetas, theta)
        it0 = it1-1
        if it1 == len(thetas):
            it1 = 0
            wt = [thetas[it1]-theta+180.0, theta-thetas[it0]]
            if r_o*2-ir0 < len(rs):
                v = arr[ir0][it0]*wr[0]*wt[0] + arr[r_o*2-ir0][it1]*wr[0]*wt[1] + arr[ir1][it0]*wr[1]*wt[0] + arr[r_o*2-ir1][it1]*wr[1]*wt[1]
            else:
                v = arr[ir0][it0]*wr[0]*wt[0] + arr[ir1][it0]*wr[1]*wt[0] # はみ出した場合、=0とする。
        elif it1 == 0:
            it0 = len(thetas)-1
            wt = [thetas[it1]-theta,theta-thetas[it0]+180.0]
            if r_o*2-ir0 < len(rs):
                v = arr[r_o*2-ir0][it0]*wr[0]*wt[0] + arr[ir0][it1]*wr[0]*wt[1] + arr[r_o*2-ir1][it0]*wr[1]*wt[0] + arr[ir1][it1]*wr[1]*wt[1]
            else:
                v = arr[ir0][it1]*wr[0]*wt[1] + arr[ir1][it1]*wr[1]*wt[1]
        else: 
            wt = [thetas[it1]-theta, theta-thetas[it0]]
            v = arr[ir0][it0]*wr[0]*wt[0] + arr[ir0][it1]*wr[0]*wt[1] + arr[ir1][it0]*wr[1]*wt[0] + arr[ir1][it1]*wr[1]*wt[1]
        return v / ((wr[0]+wr[1])*(wt[0]+wt[1]))

    # The entries of thetas are assumed to be in degrees, not radians.
    def putConvolution(self, ArrayDeriv2, image_shape, rhos, rho_o, thetas):
        # the range of rho.
        NGRID = 11
        # rhoの区間を[self.rho0_wide[0], self.rho0_wide[1]]より少し広げる。
        MinRho = self.rho0[0] - (self.rho0[0]-self.rho0_wide[0])*1.2
        MaxRho = self.rho0[1] + (self.rho0_wide[1]-self.rho0[1])*1.2
        rho_grids = np.linspace(MinRho, MaxRho, NGRID, endpoint = True)
        # print("irho:" + str(irho_begin) + ", " + str(irho_end))
        # assert irho_begin < irho_end

        # Determine the range of theta.
        LineX0=[]
        LineY0=[]
        LineX1=[]
        LineY1=[]
        getLine(image_shape, self.rho0_wide[0], self.theta0, LineX0, LineY0, self.Circle)
        getLine(image_shape, self.rho0_wide[1], self.theta0, LineX1, LineY1, self.Circle)
        if len(LineX0) < 2 or len(LineX1) < 2:
            return -1.0

        # Take the diagonals.
        image_o = [image_shape[0]//2, image_shape[1]//2]
        dtheta = 0.
        for i in range(len(LineX0)):
            p1 = [LineX0[i], LineY0[i]]
            lensq = 0.0
            index = 0
            for j in range(len(LineX1)):
                lensq2 = (LineX1[j]-p1[0])**2 + (LineY1[j]-p1[1])**2
                if lensq2 > lensq:
                    index = j
                    lensq = lensq2
            rt = getRhoThetaFromCrossings(image_o, p1, [LineX1[index], LineY1[index]])
            if rt[1] >= self.theta0 + np.pi*0.5:
                rt[1] -= np.pi
            elif rt[1] + np.pi*0.5 <= self.theta0:
                rt[1] += np.pi
            if dtheta < abs(rt[1] - self.theta0):
                dtheta = abs(rt[1] - self.theta0)
        # thetaの区間を[self.theta0 - dtheta, self.theta0 + dtheta]より少し広げる。
        MinTheta = self.theta0 - dtheta*1.2
        MaxTheta = self.theta0 + dtheta*1.2
        theta_grids = np.linspace(MinTheta, MaxTheta, NGRID, endpoint = True)

        # print("itheta:" + str(itheta_begin) + ", " + str(itheta_end))
        # assert itheta_begin < itheta_end
        ArrayMask = np.zeros((NGRID, NGRID))
        ArrayDeriv2_rescale = np.zeros((NGRID, NGRID))

        mv = 0.
        vv = 0.
        mm = 0.
        for i in range(len(rho_grids)):
            for j in range(len(theta_grids)):
                m = self.calMaskValue(image_shape, rho_grids[i], theta_grids[j])
                v = self.calRescaledValue(ArrayDeriv2, rhos, rho_o, thetas, rho_grids[i], np.rad2deg(theta_grids[j]))
                mv += m*v
                vv += v*v
                mm += m*m
                ArrayMask[i][j] = m
                ArrayDeriv2_rescale[i][j] = v
        # print(str(mv/np.sqrt(mm*vv)) +',')
        # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.5))
        # ax1.set_title('Mask')
        # ax1.set_xlabel("Projection angle (deg)")
        # ax1.set_ylabel("Projection position (pixels)")
        # ax1.imshow(ArrayMask.tolist(), cmap=plt.cm.Greys_r,
        #     extent=(MinTheta, MaxTheta, self.rho0_wide[0], self.rho0_wide[1]), aspect='auto')
        # ax2.set_title('Hough Image')
        # ax2.set_xlabel("Projection angle (deg)")
        # ax2.set_ylabel("Projection position (pixels)")
        # ax2.imshow(ArrayDeriv2_rescale.tolist(), cmap=plt.cm.Greys_r,
        #     extent=(MinTheta, MaxTheta, self.rho0_wide[0], self.rho0_wide[1]), aspect='auto')
        return mv/np.sqrt(mm*vv)

"""def demo_getLine (image_shape = (400, 400), k = 20):
    import random
    random.seed (2023)
    import itertools
    import cv2
    radius = min (image_shape) // 2
    rhos = random.choices (
        range (-radius, radius + 1), k=k)
    thetas = random.choices (
        list (np.linspace (0.0, 180.0, endpoint=False)),
        k = k)

    H, W = image_shape
    image = np.ones ((*image_shape, 3), dtype=np.uint8) * 255
    cv2.circle (image, center=(W // 2, H // 2),
                radius = radius, color = (255,255,0), thickness=3)

    #for rho, theta in itertools.product (rhos, thetas):
    for i, (rho, theta) in enumerate (zip (rhos, thetas)):
        LineX_circle, LineY_circle = [], []
        getLineCircle (image_shape, rho, theta,
                       LineX_circle, LineY_circle)
        LineX_box, LineY_box = [], []
        getLineBox (image_shape, rho, theta, LineX_box, LineY_box)
        
        plt.plot (LineX_box, LineY_box)
        plt.scatter (LineX_circle, LineY_circle,
                         #label = '{}'.format(i)
                         )
        
    plt.imshow(image)
    #plt.legend()
    plt.xlim ((0.0, W))
    plt.ylim ((0.0, H))
    plt.show()"""

#demo_getLine ((400, 400), k = 20)