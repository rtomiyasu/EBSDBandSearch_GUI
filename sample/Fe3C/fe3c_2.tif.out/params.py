PC0 = [0.4821, 0.3904, 0.5834]
# PCx,PCy,PCz (pixel)

Circle = True
# True: EBSD image is a circle, False: rectangular

RescaleParam = 400
# 画像のスケール変換

deg = 3
# 平滑化の次元

num_points = 3
# 2以上。
# 2*num_points+1点を用いて3次多項式をフィッティングする。

thred = 5.0
# ラドン変換の誤差値の計算
# (誤差値とピーク高さの比較によりバンド抽出を行う。)

MinCorrelation = 0.6
# Hough変換の2次微分のモデル値との相関は理論上[-1, 1]の範囲を取るが、
# 相関値がMinCorrelationより大きい場合にバンドとして検出する。
# Hough画像より検出したオブジェクトがバンドかどうか判定する際に使用する。

BAND_WIDTH_MIN = 0.02
BAND_WIDTH_MAX = 0.2
# バンド幅の下限と上限 (px)

dtheta = 5.0
# これより小さい角度で交わるバンドは、同一のバンドとみなす。
# Hough変換の座標(rho, theta)が近いバンドの中から最良のバンドを選定する際に使用する。



