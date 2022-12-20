
import matplotlib.pylab as plt
import numpy as np
import matplotlib
from matplotlib.pyplot import savefig
import scipy.ndimage as ndimage
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

matplotlib.rc('font', family='sans-serif')
matplotlib.rc('font', serif='Helvetica')
matplotlib.rc('text', usetex='false')
matplotlib.rcParams.update({'font.size': 6})

matplotlib.rcParams['axes.linewidth'] = 0.4

matplotlib.rcParams['xtick.major.size'] = 2.0
matplotlib.rcParams['xtick.major.width'] = 0.2
matplotlib.rcParams['xtick.minor.size'] = 1
matplotlib.rcParams['xtick.minor.width'] = 0.2

matplotlib.rcParams['ytick.major.size'] = 2.0
matplotlib.rcParams['ytick.major.width'] = 0.2
matplotlib.rcParams['ytick.minor.size'] = 1
matplotlib.rcParams['ytick.minor.width'] = 0.2

#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}",r"\renewcommand{\seriesdefault}{\bfdefault}",r"\usepackage{amsfonts}", r"\usepackage{textgreek}",r"\usepackage{textcomp}",r"\usepackage{gensymb}",r"\usepackage{fixltx2e}", r'\boldmath']
#params = {'text.usetex' : True,
#          'font.size' : 20,
#          'font.family' : 'sans-serif',
#          'text.latex.unicode': True,
#          'figure.figsize' : (8, 6), #8, 6 originally. Make it ~30% smaller so the text is ~30% bigger. (Multiply by .7)
#          'figure.autolayout' : True
#          }
#matplotlib.rcParams.update(params)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(1.18, 0.88)) # in inches
fig.subplots_adjust(left=0.18, bottom=0.22, right=0.9, top=0.96, wspace=None, hspace=None)


# CVs
#cvs = np.loadtxt('./input-cvs.txt')
#cv1=cvs[:,1]
#cv2=cvs[:,2]
#ax.scatter(cv1,cv2, s=0.5,facecolor='none', edgecolor='k')
#ax.scatter(cv1,cv2, c='k', s=0.2)


conf = "civsd-down-minus-gc-R223A-negative-500mV-yz-bias"


## PMF
a1 = np.loadtxt('./%s.txt' % (conf))

#x0data=a1[:,0]
x1data=a1[:,0]
y1data=a1[:,1]
z1data=a1[:,3]

print np.min(z1data[np.nonzero(z1data)])
print np.max(z1data[np.nonzero(z1data)])

y1list=y1data.tolist() # convert a Numpy array to a Python list
num_row=y1list.count(y1list[0])
num_column=len(y1list)/num_row

X1data=np.reshape(x1data, (num_row,num_column)) # reshape the array
Y1data=np.reshape(y1data, (num_row,num_column))
Z1data=np.reshape(z1data, (num_row,num_column))


# Now make a contour plot with the levels specified,
# and with the colormap generated automatically from a list
# of colors.
# extends = ["neither", "both", "min", "max"]
# range([start], stop[, step]) / arange() / np.linspace()


levels = np.arange(-1.4,1.5,0.2)
#levels = np.arange(-1.8,1.9,0.2)

norm=plt.Normalize(-2,2)
mycmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","white","blue"])

CS1 = ax.contourf(X1data,Y1data,Z1data,levels, alpha=1.0, cmap=mycmap, extend='neither')

#CS1 = ax.contourf(X1data,Y1data,Z1data,levels, alpha=1.0, cmap="bwr", extend='neither')

# Our data range extends outside the range of levels; make
# data below the lowest contour level yellow, and above the
# highest level cyan:

#CS1.cmap.set_under('blue')
#CS1.cmap.set_over('red')


# plot the contour levels
# We could pass in additional levels to provide extra resolution

CS11 = ax.contour(CS1, colors = 'black', linewidths = (0.1), alpha=0.5)

# Make a colorbar for the ContourSet returned by the contourf call.

#cbar = plt.colorbar(CS1)
cbar = fig.colorbar(CS1, ax=ax, ticks=[-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2], pad=0.05)

#cbar = fig.colorbar(CS1, ax=ax, ticks=[-1.0,-0.8,4,6,8,10])
#cbar = fig.colorbar(CS1, ax=ax, ticks=[0,5,10,15])

cbar.ax.set_ylabel('Electrostatic Potential (V)', fontsize=4, labelpad=2)
cbar.ax.tick_params(labelsize='2.5', pad=1)

# Add the contour line levels to the colorbar
cbar.add_lines(CS11)


## Strings

#strs0 = np.loadtxt('./string-dist-rotat-80.txt')
#strx0 = strs0[:,1]
#stry0 = strs0[:,2]
#plt.scatter(strx0, stry0, s=10, marker='s', facecolors='red', edgecolors='red')
#plt.plot(strx0, stry0, linestyle='-', color='red')


ax.set_xlabel(r'y ($\mathsf{\AA}$)', fontsize=6, labelpad=1.0)
ax.set_ylabel(r'z ($\mathsf{\AA}$)', fontsize=6, labelpad=1.0)

ax.tick_params(axis='x', labelsize= 2.5, labelcolor='black', pad=2.0)
ax.tick_params(axis='y', labelsize= 2.5, labelcolor='black', pad=2.0)


ax.set_xlim([-10, 10])
ax.set_ylim([-10, 10])
#ax.xaxis.set_ticks(np.arange(-40, 40+1, 20))
#ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
#ax.xaxis.set_minor_locator(MultipleLocator(10.0))
#ax.yaxis.set_ticks(np.arange(-40, 40+1, 20))
#ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
#ax.yaxis.set_minor_locator(MultipleLocator(10))


#plt.show()
#plt.grid(True)
savefig('%s.png' % (conf), facecolor='1.0', format='png', dpi=1200, transparent=True)
