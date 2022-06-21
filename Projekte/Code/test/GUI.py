try:
    from Tkinter import *
except ImportError:
    from tkinter import *
    
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import sys
import os
from scipy.integrate import *
from pde_int import *
from pde_functions import *
import time

def random_initialiser(xmax, ymax, NVar):
    return np.random.normal(loc=1, scale=0.1, size=xmax*ymax*NVar)

class graph():
    def __init__(self, title='Pattern Simulator', bndcondition="zeroflux", celltype="quadratic", tmax=10000):
        # simulation settings
        self.xmax=20
        self.ymax=20
        self.NVar=7
        self.ind = IJKth(1,np.arange(self.ymax),np.arange(self.xmax),self.ymax,self.NVar)
        self.shape = (self.xmax, self.ymax)
        self.k = r_[0.5982, 0.1405, 2.1971, 1.1245, 0.2916, 2.3028, 0.3466, 1.7822, 0.3976, 9.9829,\
                       1.2590, 2.6202, 1.5731, 5.2625, 4.8758, 0.3196, 0.1465, 2.1453, 0.5396, 56.0520,\
                       0.5131, 0.8396, 7.8041, 1.3647]
        self.species_labels = ('TTG1', 'GL1', 'GL3', 'TRY', 'CPC', 'AC1', 'AC2')
        self.ind=IJKth(1,np.arange(self.ymax),np.arange(self.xmax),self.ymax, self.NVar)
        self.D=couplingMatrix(self.xmax,self.ymax, bndcondition, celltype)
        self.t_span = (0, tmax)
        self.t_eval = np.linspace(self.t_span[0], self.t_span[1], 100)
        self.TD = 0
        self.CD = 0
        # GUI settings
        self.titleText = title
        self.root = Tk()
        self.root.wm_title(self.titleText)
        # w, h = self.root.winfo_screenwidth(), self.root.winfo_screenheight()
        # self.root.geometry("%dx%d+0+0" % (w, h))
        # Figure settings
        self.fig, self.ax = plt.subplots(nrows=2, ncols=4, figsize=(12,6))
        plt.tight_layout(pad=2.5)
        self.canvas = FigureCanvasTkAgg(self.fig, self.root)
        self.canvas.get_tk_widget().grid(row=0, column=1)
        # Labels and entries
        self.parlabels = ['k'+str(x) for x in range(1, len(self.k)+1)]
        self.labels = []
        self.entries = []
        group = LabelFrame(self.root, text = "Parameters", padx = 5, pady = 5)
        group.grid(row=0, column=0)
        # First column of k1 to k12
        for i, l in enumerate(self.parlabels):
            label = Label(group, text=l)
            entry = Entry(group, width = 6)
            if i < len(self.parlabels)/2:
                label.grid(row=i, column=0)
                entry.grid(row=i, column=1)
            else:
                label.grid(row=i-len(self.parlabels)//2, column=2)
                entry.grid(row=i-len(self.parlabels)//2, column=3)
            entry.insert(10, self.k[i])
            self.labels.append(label)
            self.entries.append(entry)
        infoframe = LabelFrame(self.root, text = "Pattern data", padx = 5, pady = 5)
        infoframe.grid(row = 1, column = 0)
        labelTD1 = Label(infoframe, text = 'Trichome density: ')
        labelTD1.grid(row = 0, column = 0)
        self.labelTD2 = Label(infoframe, text = '{:.2f} %'.format(self.TD))
        self.labelTD2.grid(row = 0, column = 1)
        labelCD1 = Label(infoframe, text = 'Cluster density: ')
        labelCD1.grid(row=1, column = 0)
        self.labelCD2 = Label(infoframe, text = '{:.2f} %'.format(self.CD))
        self.labelCD2.grid(row =1 , column = 1)

        # Buttons
        self.buttonFrame = Frame()
        self.buttonFrame.grid(row=1,column=1)
        run = Button(self.buttonFrame, text = "Run", command = self.run)
        run.pack(side="left",fill='x')
        reinit = Button(self.buttonFrame, text="Re-initialize", command = self.init)   
        reinit.pack(side="left",fill='x')
        reset = Button(self.buttonFrame, text="Reset parameters", command = self.reset)
        reset.pack(side="left", fill='x')
        quit = Button(self.buttonFrame, text = "Quit", command = self.root.quit)
        quit.pack(side="left",fill='x')


    def draw(self,ax,xmax,ymax,y,ind):
        ax = ax.flatten()
        X,Y=meshgrid(np.arange(xmax),np.arange(ymax))
        for i in range(self.NVar):
            ax[i].clear()
            # Major ticks
            ax[i].set_xticks(np.arange(-0.55, self.xmax, 1.0));
            ax[i].set_yticks(np.arange(-0.6, self.ymax, 1.0));

            ax[i].set_xticklabels([])
            ax[i].set_yticklabels([])
            # Gridlines based on minor ticks
            ax[i].grid(which='major', color='black', linestyle='-', linewidth=1)
            ax[i].set_title(self.species_labels[i])
            # ax[i].get_xaxis().set_visible(False)
            # ax[i].get_yaxis().set_visible(False)
            # ax[i].contourf(X,Y,y[ind+i].reshape(shape(X)),30,cmap='YlGn')
            ax[i].imshow(y[ind+i].reshape(shape(X)), cmap='YlGn')
        ax[self.NVar].clear()
        # Major ticks
        ax[self.NVar].set_xticks(np.arange(-0.55, self.xmax, 1.0));
        ax[self.NVar].set_yticks(np.arange(-0.6, self.ymax, 1.0));

        ax[self.NVar].set_xticklabels([])
        ax[self.NVar].set_yticklabels([])
        # Gridlines based on minor ticks
        ax[self.NVar].grid(which='major', color='black', linestyle='-', linewidth=1)
        ax[self.NVar].set_title('AC1+AC2')
        # ax[self.NVar].get_xaxis().set_visible(False)
        # ax[self.NVar].get_yaxis().set_visible(False)
        # ax[self.NVar].contourf(X,Y,(y[ind+5] + y[ind+6]).reshape(shape(X)),30,cmap='YlGn')
        ax[self.NVar].imshow((y[ind+5] + y[ind+6]).reshape(shape(X)), cmap='YlGn')
        self.canvas.draw()

    def pattern_quant(self, y, ind):
        # Trichome density
        ACsum = y[ind+5] + y[ind+6];
        ACmax = amax(ACsum)
        T = np.nonzero(ACsum > 0.5*ACmax)
        TD = (T[0].size/400) * 100
        # Cluster density
        j,i = np.nonzero(self.D)
        inCluster = np.zeros(self.xmax*self.ymax)
        for Tidx in T[0]:
            nidx = i[j==Tidx]
            nidx = nidx[nidx != Tidx]
            ismember = [nb in T[0] for nb in nidx]
            if any(ismember):
                inCluster[Tidx] = 1
        CD = (np.sum(inCluster)/T[0].size)*100
        return TD, CD



    def get_parameters(self):
        self.pars = []
        for entry in self.entries:
            self.pars.append(float(entry.get()))
        return self.pars

    def start(self):
        self.init()
        self.root.mainloop()

    def run(self):
        t0=time.time()
        x = self.root.winfo_x()
        y = self.root.winfo_y()
        window = Toplevel(self.root)
        window.title("Message")
        label = Label(window, text = "Calculating... Please wait.", padx=100, pady=100)
        label.pack()
        dx = self.root.winfo_width()//2
        dy = self.root.winfo_height()//2

        window.geometry("%dx%d+%d+%d" % (200, 100, x+dx, y+dy))
        window.update()
        print("Integration started")
        self.k = self.get_parameters()
        sol = solve_ivp(lambda t, y: full_model(t, y, self.D, self.ind, self.k), self.t_span,\
                        self.y0, method='Radau',jac_sparsity=jpat(self.D,self.ind), vectorized=True, t_eval=self.t_eval)
        print("Integration done in %(sec)d seconds" %{"sec":(time.time()-t0)})
        window.destroy()
        yss = sol.y[:, -1]
        self.draw(self.ax, self.xmax, self.ymax, yss, self.ind)
        self.TD, self.CD = self.pattern_quant(yss, self.ind)
        self.labelTD2['text'] = '{:.2f} %'.format(self.TD)
        self.labelCD2['text'] = '{:.2f} %'.format(self.CD)    


    def init(self):
        self.k = self.get_parameters()
        self.t = 0
        self.y0 = random_initialiser(self.xmax, self.ymax, self.NVar)
        self.draw(self.ax, self.xmax, self.ymax, self.y0, self.ind)
        self.TD = 0
        self.CD = 0
        self.labelTD2['text'] = '{:.2f} %'.format(self.TD)
        self.labelCD2['text'] = '{:.2f} %'.format(self.CD) 

    def reset(self):
        self.k=r_[0.5982, 0.1405, 2.1971, 1.1245, 0.2916, 2.3028, 0.3466, 1.7822, 0.3976, 9.9829,\
                       1.2590, 2.6202, 1.5731, 5.2625, 4.8758, 0.3196, 0.1465, 2.1453, 0.5396, 56.0520,\
                       0.5131, 0.8396, 7.8041, 1.3647]
        for i, entry in enumerate(self.entries):
            entry.delete(0, 'end')
            entry.insert(10, self.k[i])

G = graph()
G.start()

