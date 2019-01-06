from numpy import sin,cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
import time

class DoublePendulum:

    def __init__(self,init_state,
                    l1,l2,
                    m1, m2,
                    g, dens,
                    origin,water):

        self.l1=l1
        self.l2=l2

        if water=='Yes':
            m1*=(1+0.5*(1/dens))
            m2*=(1+0.5*(1/dens))
            g*=(1-1/dens)
            #const self.g
        self.m1=m1
        self.m2=m2
        self.g=g
        self.dens=dens
        print self.g
        self.init_state=np.asarray(init_state,dtype='float')
        self.params=(self.l1,self.l2,self.m1,self.m2,g)
        self.origin=origin
        self.time_elapsed=0
        self.state=self.init_state*np.pi/180

    def position(self):
        (self.l1,self.l2,self.m1,self.m2,self.g)=self.params
        x=np.cumsum([self.origin[0],
                    self.l1*sin(self.state[0]),
                    self.l2*sin(self.state[2])])

        y=np.cumsum([self.origin[1],
                    -self.l1*cos(self.state[0]),
                    -self.l2*cos(self.state[2])])
        return (x,y)


    def energy(self):
        (self.l1,self.l2,self.m1,self.m2,self.g)=self.params
        x=np.cumsum([
                    self.l1*sin(self.state[0]),
                    self.l2*sin(self.state[2])])
        y=np.cumsum([
                    -self.l1*cos(self.state[0]),
                    -self.l2*cos(self.state[2])])
        vx=np.cumsum([self.l1*self.state[1]*cos(self.state[0]),
                    self.l2*self.state[3]*cos(self.state[2])])

        vy=np.cumsum([self.l1*self.state[1]*sin(self.state[0]),
                    self.l2*self.state[3]*sin(self.state[2])])

        u=self.g*(self.m1*y[0]+self.m2*y[1])
        k=0.5*(self.m1*(vx[0]**2+vy[0]**2)+self.m2*(vx[1]**2+vy[1]**2))

        return (u,k)

    def dstate_dt(self,state,t):
        (self.l1,self.l2,self.m1,self.m2,self.g)=self.params
        dxdt=np.asarray([0,0,0,0],float)
        dxdt[0]=state[1]
        dxdt[2]=state[3]
        cos_delta=cos(state[0]-state[2])
        sin_delta=sin(state[0]-state[2])

        num1=-self.g*(2*self.m1+self.m2)*sin(state[0])
        num2=-self.m2*self.g*sin(state[0]-state[2]*2)
        num3=-2*sin_delta*self.m2
        num4=state[3]*state[3]*self.l2+state[1]*state[1]*self.l1*cos_delta
        den1=self.l1*(2*self.m1+self.m2-self.m2*cos(2*state[0]-2*state[2]))

        dxdt[1]=(num1+num2+num3*num4)/den1
        #print self.g,dxdt[1]
        num1=2*sin(state[0]-state[2])
        num2=state[1]*state[1]*self.l1*(self.m1+self.m2)
        num3=self.g*(self.m1+self.m2)*cos(state[0])
        num4=state[3]*state[3]*self.l2*self.m2*cos_delta
        den2=self.l2*(2*self.m1+self.m2-self.m2*cos(2*state[0]-2*state[2]))
        dxdt[3]=(num1*(num2+num3+num4))/den2

        return dxdt

    def step(self,dt,ax):
        self.state=integrate.odeint(self.dstate_dt,self.state,[0,dt])[1]
        self.time_elapsed+=dt
        self.l1,self.l2=self.params[0],self.params[1]
        x=np.cumsum([
                    self.l1*sin(self.state[0]),
                    self.l2*sin(self.state[2])])
        y=np.cumsum([
                    -self.l1*cos(self.state[0]),
                    -self.l2*cos(self.state[2])])
        plt.sca(ax)
        plt.scatter(x[0],y[0],s=5,c='red')
        plt.scatter(x[1],y[1],s=5,c='pink')


#pendulum=DoublePendulum([90.,0.0,90.,0.0],3,5,2,10,7.924,(0,0),'No')
pendulum=DoublePendulum([90.,0.0,90.,0.0],3,5,2,10,19.8,7.86,(0,0),'Yes')
pendulum2=DoublePendulum([90.,0.0,90.,0.0],3,5,2,10,19.8,7.86,(0,0),'No')
dt=1./25

fig=plt.figure(figsize=(80,80),facecolor= u'#191970')
ax=fig.add_subplot(2,2,2,aspect='equal',autoscale_on=False,xlim=(-10,10),ylim=(-10,10))

img=plt.imread("water2.jpg")
#plt.xlabel('x',color='w')
# plt.ylabel('y')
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
#plt.xlabel('x1',color='w')
plt.rc_context({'axes.edgecolor':'orange', 'xtick.color':'w', 'ytick.color':'w'})
ax.tick_params(axis='x', colors='w')
ax.tick_params(axis='y', colors='w')

time_text=ax.text(0.02,0.95,'',transform=ax.transAxes,bbox=dict(facecolor='w', alpha=0.5),fontsize='12')
time_text2=ax2.text(0.02,0.95,'',transform=ax.transAxes,bbox=dict(facecolor='w', alpha=0.5),fontsize='12')
l1_text=ax.text(0.02,0.95,'',transform=ax.transAxes,bbox=dict(facecolor='w', alpha=0.5),fontsize='6')
l2_text=ax2.text(0.02,0.95,'',transform=ax.transAxes,bbox=dict(facecolor='w', alpha=0.5),fontsize='6')
m1_text=ax.text(0.02,0.95,'',transform=ax.transAxes)
m2_text=ax2.text(0.02,0.95,'',transform=ax.transAxes)
energy_text=ax.text(0.02,0.85,'',transform=ax.transAxes,bbox=dict(facecolor='w', alpha=0.5),fontsize='12')
energy_text2=ax2.text(0.02,0.85,'',transform=ax.transAxes,bbox=dict(facecolor='w', alpha=0.5),fontsize='12')

def init():
    line.set_data([],[])
    line2.set_data([],[])
    time_text.set_text('')
    energy_text.set_text('')
    time_text2.set_text('')
    energy_text2.set_text('')
    return line,line2,l1_text,l2_text,m1_text,m2_text,energy_text,energy_text2,time_text2


def animate(i):
    global pendulum,pendulum2,dt
    pendulum2.step(dt,ax2)
    pendulum.step(dt,ax)
    line.set_data(*pendulum.position())
    line2.set_data(*pendulum2.position())
    plt.sca(ax)
    time_text.set_text('Time=%.1fs' % pendulum.time_elapsed)
    energy_text.set_text('\nKinetic energy= %.3f J' % (pendulum.energy()[1]))
    plt.sca(ax)
    time_text2.set_text('Time=%.1fs' % pendulum2.time_elapsed)
    energy_text2.set_text('\nKinetic energy= %.3f J' % pendulum2.energy()[1])
    return line,line2,time_text,energy_text,time_text2,energy_text2


plt.title('PENDULUM IN AIR\n',color='w',fontsize=17)
plt.sca(ax)
plt.title('PENDULUM IN WATER\n',color='w',fontsize='17')
textstr = '\n\nLength of pendulum1=%.2f m \n\nLength of pendulum2=%.2f m\n\nMass of pendulum1=%.2f kg \n\nMass of pendulum2=%.2f kg\n\nAcceleration due to gravity=%.2f N/kg\n\nDensity of bob=%.2f g/cm**3\n\nDensity of water=1.00g/cm**3'%(pendulum.l1,pendulum.l2,pendulum2.m1,pendulum2.m2,pendulum2.g,pendulum.dens)
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
