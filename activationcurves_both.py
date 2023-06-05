import matplotlib.pyplot as plt
import numpy as np

def mInf_graph(v,x,y):
    N = len(v)
    m0 = 0 
    mInf   = np.zeros(N)
    mTau   = np.zeros(N)
    m      = np.zeros(N)
    for i in range(N):
        xi = x[i]
        if abs(xi/y)<1e-6:
            mAlpha = 0.055*y*(1-xi)/(2.*y) 
        else:
            mAlpha = 0.055*xi/(np.exp(xi/y)-1)
        mBeta = 0.94*np.exp((-75-v[i])/17.)
        mInf[i] = mAlpha/(mAlpha+mBeta)
        mTau[i] = 1./(mAlpha+mBeta)
        if i==0:
            m[i]    = mInf[i]-(mInf[i]-m0)*np.exp(-300/mTau[i])
        else:
            m[i]    = mInf[i]-(mInf[i]-m[i-1])*np.exp(-300/mTau[i])
    return m, mInf, mTau

def hInf_graph(v):
    N = len(v)
    h0 = 0 
    hInf   = np.zeros(N)
    hTau   = np.zeros(N)
    h      = np.zeros(N)
    for i in range(N):
        vi = v[i]
        hAlpha = 0.000457*np.exp(-(13+vi)/50.)
        hBeta = 0.0065/(np.exp(-(15+vi)/28.)+1)
        hInf[i] = hAlpha/(hAlpha+hBeta)
        hTau[i] = 1./(hAlpha+hBeta)
        if i==0:
            h[i]    = hInf[i]-(hInf[i]-h0)*np.exp(-300/hTau[i])
        else:
            h[i]    = hInf[i]-(hInf[i]-h[i-1])*np.exp(-300/hTau[i])
    return h, hInf, hTau

def zInf_graph(cai):
    N = len(cai)
    zInf = np.zeros(N)
    zInf = 1/(1+(0.00043/cai)**4.8)
    return zInf

cais = np.logspace(-5, -2, num=100, endpoint=True, base=10.0, dtype=None, axis=0)
zInf_standard = zInf_graph(cais)


caihalf_SK = 0
passedit_SK = False
for i in range(len(zInf_standard)):
    if zInf_standard[i]>0.5 and passedit_SK==False:
        caihalf_SK = cais[i]
        passedit_SK=True
        break

vs = np.linspace(-150,150,1000)

xstandard = -27-vs
ystandard = 3.8
m_standard, mInf_standard, mTau_standard = mInf_graph(vs,xstandard,ystandard)
h_standard, hInf_standard, hTau_standard = hInf_graph(vs)

Vhalf = 0
passedit = False
for i in range(len(mInf_standard)):
    if mInf_standard[i]>0.5 and passedit==False:
        Vhalf = vs[i]
        passedit=True
        break

print('max(mTau_standard):',max(mTau_standard),'; min(mTau_standard):',min(mTau_standard))
print('max(hTau_standard):',max(hTau_standard),'; min(hTau_standard):',min(hTau_standard))

print('caihalf_SK:',caihalf_SK)
plt.figure(figsize=(5,2.7),dpi=300)
plt.plot(cais,zInf_standard,color='k')
plt.axvline(x=1e-4,color='grey',linestyle=':',linewidth=0.75)
plt.axvline(x=caihalf_SK,color='k',linestyle='--',linewidth=0.75)
plt.xscale("log")
plt.xlabel(r'[Ca$^{2+}$]$_\mathregular{in}$ (mM)',fontsize=12)
plt.ylabel(r'Activation',fontsize=12)
plt.title(r'$z_\infty$',fontsize=18)
plt.tight_layout()
plt.savefig('zInf_SK.png')
plt.show()

