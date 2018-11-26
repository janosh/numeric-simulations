import matplotlib.pyplot as plt
import numpy as np

def frame(phi1, phi2, phi1_arr, phi2_arr, t, num):
    #length of the pendula
    L1 = 2.
    L2 = 1.
    L = L1 + L2

    #create a new figure
    fig = plt.figure(figsize=[8,8], facecolor='white')
    axes = fig.gca()

    #proportional x and y axis
    axes.axis("equal")
    #no coordinate labels
    axes.set_axis_off()
    #set the plot range in x and y direction
    axes.set_ylim([-(L+0.2),(L+0.2)])
    axes.set_xlim([-(L+0.2),(L+0.2)])


    #calculate position of pendula
    pos1 = np.zeros(2)
    pos2 = np.zeros(2)

    pos1[0] = L1 * np.sin(phi1)
    pos1[1] = - L1 * np.cos(phi1)
    pos2[0] = pos1[0] + L2 * np.sin(phi2)
    pos2[1] = pos1[1] - L2 * np.cos(phi2)

    #plot the two pendula
    c0 = plt.Circle([0.,0.],.07, color='k', zorder=100)
    c1 = plt.Circle(pos1,.05, color='r', zorder=100)
    c2 = plt.Circle(pos2,.05, color='b', zorder=100)
    axes.add_artist(c0)
    axes.add_artist(c1)
    axes.add_artist(c2)

    #plot the two limbs
    axes.plot([0., pos1[0]], [0., pos1[1]], color='k', linewidth=2, zorder=50)
    axes.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]], color='k', linewidth=2, zorder=50)

    #add the trajectory in the background
    N = len(phi1_arr)
    pos1_arr = np.zeros((N,2))
    pos2_arr = np.zeros((N,2))

    pos1_arr[:,0] = L1 * np.sin(phi1_arr)
    pos1_arr[:,1] = - L1 * np.cos(phi1_arr)
    pos2_arr[:,0] = pos1_arr[:,0] + L2 * np.sin(phi2_arr)
    pos2_arr[:,1] = pos1_arr[:,1] - L2 * np.cos(phi2_arr)

    axes.plot(pos1_arr[:,0], pos1_arr[:,1], color='r', linewidth = 0.5, zorder=0)
    axes.plot(pos2_arr[:,0], pos2_arr[:,1], color='b', linewidth = 0.5, zorder=0)

    #add a time label
    axes.text(0.7,0.85, "T = %03.2f"%t, transform = fig.transFigure, fontsize=16, horizontalalignment='left', verticalalignment='center')

    #save the figure
    fig.savefig("Frame Buffer/frame_%03d.png"%num, dpi=80)

    plt.close(fig)


if __name__ == "__main__":
    data_e = np.fromfile('movie.dat', dtype=float, count=-1, sep='')
    phi1 = data_e[::2]
    phi2 = data_e[1::2]

    for i in np.arange(2000):
        frame(phi1[i], phi2[i], phi1[:i+1], phi2[:i+1], i*0.2, i)
