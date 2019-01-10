from numpy import sin,cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
import time

class Pendulum:

    def __init__(self,init_condition,
                    l1,l2,
                    m1, m2,
                    g, dens,water):
        self.l1=l1
        self.l2=l2
        if water=='Yes':
            m1*=(1+0.5*(1/dens))
            m2*=(1+0.5*(1/dens))
            g*=(1-1/dens)
        self.m1=m1
        self.m2=m2
        self.g=g
        self.dens=dens
        self.condition=np.asarray(init_condition,dtype='float')
        self.params=(self.l1,self.l2,self.m1,self.m2,g)
        self.time_elapsed=0
        self.initial_coordinate=(0.0,8.0)
        self.condition*=np.pi/180

    def find_posn(self):
        (self.l1,self.l2,self.m1,self.m2,self.g)=self.params
        x=np.cumsum([0,
                    self.l1*sin(self.condition[0]),
                    self.l2*sin(self.condition[2])])

        y=np.cumsum([0,
                    -self.l1*cos(self.condition[0]),
                    -self.l2*cos(self.condition[2])])

        return (x,y)


    def find_energy(self):
        (self.l1,self.l2,self.m1,self.m2,self.g)=self.params
        x=np.cumsum([
                    self.l1*sin(self.condition[0]),
                    self.l2*sin(self.condition[2])])
        y=np.cumsum([
                    -self.l1*cos(self.condition[0]),
                    -self.l2*cos(self.condition[2])])
        vx=np.cumsum([self.l1*self.condition[1]*cos(self.condition[0]),
                    self.l2*self.condition[3]*cos(self.condition[2])])

        vy=np.cumsum([self.l1*self.condition[1]*sin(self.condition[0]),
                    self.l2*self.condition[3]*sin(self.condition[2])])

        total_energy=self.g*(self.m1*self.initial_coordinate[0]+self.m2*self.initial_coordinate[1])

        new=integrate.odeint(self.differentiate,self.condition,[0,dt])[1]
        i1=self.m1*self.l1*self.l1
        i2=self.m2*self.l2*self.l2
        q=0.5*(i1*new[0]*new[0] + i2*new[2]*new[2])
        k=0.5*(self.m1*(vx[0]**2+vy[0]**2)+self.m2*(vx[1]**2+vy[1]**2))
        k+=q
        u=total_energy-k

        return (u,k)

    def differentiate(self,condition,t):
        (self.l1,self.l2,self.m1,self.m2,self.g)=self.params
        dxdt=np.asarray([0,0,0,0],float)
        dxdt[0]=condition[1]
        dxdt[2]=condition[3]

        a=-self.g*(2*self.m1+self.m2)*sin(condition[0])
        b=-self.m2*self.g*sin(condition[0]-condition[2]*2)
        c=-2*sin(condition[0]-condition[2])*self.m2
        d=condition[3]*condition[3]*self.l2+condition[1]*condition[1]*self.l1*cos(condition[0]-condition[2])
        e=self.l1*(2*self.m1+self.m2-self.m2*cos(2*condition[0]-2*condition[2]))

        dxdt[1]=(a+b+c*d)/e
        a=2*sin(condition[0]-condition[2])
        b=condition[1]*condition[1]*self.l1*(self.m1+self.m2)
        c=self.g*(self.m1+self.m2)*cos(condition[0])
        d=condition[3]*condition[3]*self.l2*self.m2*cos(condition[0]-condition[2])
        e=self.l2*(2*self.m1+self.m2-self.m2*cos(2*condition[0]-2*condition[2]))
        dxdt[3]=(a*(b+c+d))/e

        return dxdt

    def increment_time(self,dt,ax):
        self.condition=integrate.odeint(self.differentiate,self.condition,[0,dt])[1]
        self.time_elapsed+=dt
        self.l1,self.l2=self.params[0],self.params[1]
        x=np.cumsum([
                    self.l1*sin(self.condition[0]),
                    self.l2*sin(self.condition[2])])
        y=np.cumsum([
                    -self.l1*cos(self.condition[0]),
                    -self.l2*cos(self.condition[2])])
        plt.sca(ax)
        plt.scatter(x[0],y[0],s=5,c='red')
        plt.scatter(x[1],y[1],s=5,c='pink')

p1=Pendulum([90.,0.0,90.,0.0],3,5,1,5,19.8,7.86,'Yes')     #CHANGE THESE VALUES TO SEE VARIATION IN BEHAVIOUR OF PENDULUM
p2=Pendulum([90.,0.0,90.,0.0],3,5,1,5,19.8,7.86,'No')
dt=1./25

fig=plt.figure(figsize=(80,80),facecolor= u'#191970')
ax=fig.add_subplot(2,2,2,aspect='equal',autoscale_on=False,xlim=(-10,10),ylim=(-10,10))
img=plt.imread("water2.jpg")
ax.imshow(img,extent=[-10,10,-10,25])
ax.grid()
plt.rc_context({'axes.edgecolor':'orange', 'xtick.color':'w', 'ytick.color':'w'})

line, = ax.plot([],'o-',color='black',markersize='15',lw=2)

ax2=fig.add_subplot(2,2,1,aspect='equal',autoscale_on=False,xlim=(-10,10),ylim=(-10,10))
ax2.grid()
plt.rc_context({'axes.edgecolor':'orange', 'xtick.color':'w', 'ytick.color':'w'})
ax2.set_facecolor('black')
line2, = ax2.plot([],'bo-',markersize='15',lw=2)
plt.sca(ax2)
plt.rc_context({'axes.edgecolor':'orange', 'xtick.color':'w', 'ytick.color':'w'})
ax.tick_params(axis='x', colors='w')
ax.tick_params(axis='y', colors='w')

time_data=ax2.text(0.72,0.95,'',transform=ax.transAxes,bbox=dict(facecolor='w', alpha=0.5),fontsize='12')
energy_data=ax2.text(0.02,0.8,'',transform=ax.transAxes,fontsize='12')
energy_data2=ax2.text(0.02,0.8,'',transform=ax2.transAxes,color='w',fontsize='12')

def init():
    line.set_data([],[])
    line2.set_data([],[])
    time_data.set_text('')
    energy_data.set_text('')
    energy_data2.set_text('')
    return line,line2,
    energy_data,energy_data2,time_data


def animate(i):
    global p1,p2,dt
    p2.increment_time(dt,ax2)
    p1.increment_time(dt,ax)
    line.set_data(*p1.find_posn())
    line2.set_data(*p2.find_posn())

    time_data.set_text('Time=%.1fs' % p1.time_elapsed)
    energy_data.set_text('\nKinetic energy= %.3f J\nPotential energy= %3f J\n\nTotal energy= %.3f J'%(p1.find_energy()[1], p1.find_energy()[0],p1.find_energy()[1]+ p1.find_energy()[0]))
    energy_data2.set_text('\n\nKinetic energy= %.3f J\nPotential energy= %3f J\n\nTotal energy= %.3f J'%(p2.find_energy()[1], p2.find_energy()[0],p2.find_energy()[1]+ p2.find_energy()[0]))
    return line,line2,
    time_data,
    energy_data,energy_data2

plt.title('PENDULUM IN AIR\n',color='w',fontsize=17)
plt.sca(ax)
plt.title('PENDULUM IN WATER\n',color='w',fontsize='17')
textstr = '\n\nLength of pendulum1=%.2f m \n\nLength of pendulum2=%.2f m\n\nMass of pendulum1=%.2f kg \n\nMass of pendulum2=%.2f kg\n\nAcceleration due to gravity=%.2f N/kg\n\nDensity of bob=%.2f g/cm**3\n\nDensity of water=1.00g/cm**3'%(p1.l1,p1.l2,p2.m1,p2.m2,p2.g,p1.dens)
plt.gcf().text(0.03, 0.5, textstr, fontsize=13,color='w')
plt.gcf().canvas.set_window_title('Double Pendulum Simulation')
plt.subplots_adjust(left=.25,top=.95,bottom=-.9)

from time import time
t0=time()
animate(0)
t1=time()
interval=1000*dt-(t1-t0)

ani=animation.FuncAnimation(fig,animate,frames=150,interval=interval,blit=False,init_func=init)
plt.show()
