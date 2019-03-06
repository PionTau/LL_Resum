import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

# Get the values for Cuhre
uM_val = np.loadtxt("PlotResumFinal_13_LL.dat", unpack = True)[1]
result = np.loadtxt("PlotResumFinal_13_LL.dat", unpack = True)[2]
errorR = np.loadtxt("PlotResumFinal_13_LL.dat", unpack = True)[3]

# Ge the values for Vegas
uM_valCSS = np.loadtxt("data_ggHiggsCuhre.text", unpack = True)[0]
resultCSS = np.loadtxt("data_ggHiggsCuhre.text", unpack = True)[1]
errorRCSS = np.loadtxt("data_ggHiggsCuhre.text", unpack = True)[2]

# Plot the figure
fig = plt.figure()

plt.plot(uM_val, result, 'b',label="LLgg Mu", linewidth=1.5)
plt.plot(uM_valCSS, resultCSS, 'r--',label="LLgg Me", linewidth=1.5)
plt.title(r'Small $p_T$ Resummation. ($\sqrt{s}=13 TeV, m_H=125 GeV$)')
# plt.ylim(bottom=-0.026)
plt.xlabel('pt')
plt.ylabel(r'$\frac{d\sigma}{dp_T}(pb/GeV)$')
# plt.grid(True)
# plt.fill_between(uM_valCuhre, resultCuhre-errorRCuhre, resultCuhre+errorRCuhre, color='gray', alpha=0.2)
plt.legend()

# Adjust layout
fig.tight_layout()

fig.savefig('PlotResumFinal_13_LL.png', dpi=100)