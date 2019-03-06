import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

# Get the values for Cuhre
uM_val = np.loadtxt("PlotResumFinal_13_LL_tiny.dat", unpack = True)[1]
result1 = np.loadtxt("PlotResumFinal_13_LL_tiny.dat", unpack = True)[2]
result2 = np.loadtxt("PlotResumFinal_13_LL_tiny.dat", unpack = True)[6]
result3 = np.loadtxt("PlotResumFinal_13_LL_tiny.dat", unpack = True)[8]

# # Ge the values for Vegas
# uM_valCSS = np.loadtxt("data_ggHiggsCuhre.text", unpack = True)[0]
# resultCSS = np.loadtxt("data_ggHiggsCuhre.text", unpack = True)[1]
# errorRCSS = np.loadtxt("data_ggHiggsCuhre.text", unpack = True)[2]

# Plot the figure
fig = plt.figure()

plt.plot(uM_val, result1, 'b',label="Small pt", linewidth=1.5)
plt.plot(uM_val, result2, 'r--',label="Cons pt", linewidth=1.5)
plt.plot(uM_val, result3, 'black',label="CSS", linewidth=1.5)
plt.title(r'Small $p_T$ Resummation. ($\sqrt{s}=13 TeV, m_H=125 GeV$)')
plt.ylim(top=1)
plt.xlabel('pt')
plt.ylabel(r'$\frac{d\sigma}{dp_T}(pb/GeV)$')
# plt.grid(True)
# plt.fill_between(uM_valCuhre, resultCuhre-errorRCuhre, resultCuhre+errorRCuhre, color='gray', alpha=0.2)
plt.legend()

# Adjust layout
fig.tight_layout()

fig.savefig('PlotResumFinal_13_LL2.png', dpi=100)