import rebound
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns
import csv
import matplotlib
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation
import glob

from mpl_toolkits.mplot3d import Axes3D
from IPython.display import clear_output
from tqdm import tqdm




class Particle:
    """
    Particle Class.

    Parameters
    ----------
    m: scalar
        Mass of the particle in solar masses.
    x,y,z: scalars
        Distances of the particle from the origin in AU.
    vx,vy,vz: scalars
        Initial velocity components of the particle in AU/yr.
    """
    def __init__(self,m,x=0,y=0,z=0,vx=0,vy=0,vz=0):
        self.m = m
        self.x = x
        self.y = y
        self.z = z
        self.a = np.linalg.norm([x, y, z])
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.vi = np.linalg.norm([vx, vy, vz])


class StarSystem:
    G = 39.476926421373

    """
    Star System Class.

    Parameters
    ----------
    stars: list
        List of star objects.
    ast: particle object
        Particle object.
    n_ast: int
        Number of asteroids.
    seed: int
        Seed with which the asteroids initial position
        angles will be randomly generated.
    """

    def __init__(self,stars,ast,n_ast,seed=0):
        self.stars = stars
        self.ast = ast
        self.n_ast = n_ast
        self.seed = seed


    def spawnAsteroid(self):
        global G
        spherical_to_cartesian = lambda rho,phi,theta : [rho*np.sin(phi)*np.cos(theta),
                                                          rho*np.cos(phi)*np.sin(theta),
                                                          rho*np.cos(phi)]
        random_angles = lambda : [2*np.pi*np.random.random_sample(),np.pi*np.random.random_sample()]
        [theta,phi] = random_angles()
        r = spherical_to_cartesian(self.ast.a,phi,theta)
        M = 0
        for star in self.stars:
            M += star.m
        rho = np.random.normal(0,(1/3)*np.sqrt(G*M/self.ast.a),1)
#         rho = np.sqrt(G*M/self.ast.a) * np.random.random_sample()
        while (rho < 0 or rho > np.sqrt(G*M/self.ast.a)):
#             rho = np.sqrt(G*M/self.ast.a) * np.random.random_sample()
            rho = np.random.normal(0,(1/3)*np.sqrt(G*M/self.ast.a),1)
        print(type(r))
        print(type(-rho/np.linalg.norm(r)))
        v_r = (-rho/np.linalg.norm(r))*r
        k = np.sqrt(G*M/self.ast.a - np.linalg.norm(v_r)**2)
        m = 2*np.random.random_sample() - 1
        n = 2*np.random.random_sample() - 1
        v_t = m*np.array([1,0,-r[0]/r[2]]) + n*np.array([0,1,-r[1]/r[2]])
        v_t = (k/np.linalg.norm(v_t))*v_t
        v = v_r + v_t
        return r,v


    def initiateSim(self):
        """
        Creates simulation in Rebound.

        Returns
        ----------
        sim: Rebound simulation object
            Simulation object with parameters of the BinaryStarClass.
        """
        sim = rebound.Simulation()
        sim.units = ("yr","AU","Msun")
        sim.G = 39.476926421373

        for star in self.stars:
            sim.add(m=star.m,
                    x=star.x,y=star.y,z=star.z,
                    vx=star.vx,vy=star.vy,vz=star.vz)

        for i in range(self.n_ast):
            ast = self.spawnAsteroid()
            sim.add(m=self.ast.m,x=ast[0][0],y=ast[0][1],z=ast[0][2],vx=ast[1][0],vy=ast[1][1],vz=ast[1][2])
        return sim


    def integrateSim(self,t_start,t_stop,iterations,sim,dt="reverselog",logbase=np.e):
        """
        Integrates simulation for a set time and iteration count.

        Parameters
        ----------
        t: scalar
            End time of the integration in years
        iterations: int
            Number of times the integration is performed
        sim: Rebound.Simulation() object
            Rebound simulation
        dt: String
            Specifies whether to integrate linearly or logarithmically
        logbase: scalar
            Specifies the base of the log, should the user choose logarithmic integration

        Returns
        ----------
        data: list
            List of dataframes, each containing the t,x,y,z,d,vx,vy,vz,v
            of each star and asteroid for every integration iteration.
        """
        global G
        mass = 0
        for i in range(len(self.stars)):
            mass += self.stars[i].m
        for i in range(self.n_ast):
            mass += self.ast.m
        print(mass)
        data = [np.zeros((iterations,11)) for i in range(sim.N)]

        if dt == "reverselog":
            timespace = t_start+t_stop-np.flipud(np.logspace(start=0,
                                                stop=np.log(t_stop)/np.log(logbase),
                                                num=iterations,
                                                base=logbase))
        elif dt == "linear":
            timespace = np.linspace(start=t_start,stop=t_stop,num=iterations)
        else:
            dt = input("Enter \"reverselog\" or \"linear\": ")
            return self.integrateSim(t_start,t_stop,iterations,sim,dt)

        i = 0
        for time in tqdm(timespace):
            clear_output(wait=True)
            sim.integrate(time)
            for j in range(sim.N):
                data[j][i,0:1] = time
                data[j][i,1:4] = sim.particles[j].xyz
                data[j][i,4:5] = np.linalg.norm(data[j][i,1:4])
                data[j][i,5:8] = np.array(sim.particles[j].vxyz)#*4.74372
                data[j][i,8:9] = np.linalg.norm(data[j][i,5:8])



                d_1 = np.sqrt((data[j][i,1:2] - data[0][i,1:2])**2 + (data[j][i,2:3] - data[0][i,2:3])**2 + (data[j][i,3:4] - data[0][i,4:5])**2)
                d_2 = np.sqrt((data[j][i,1:2] - data[1][i,1:2])**2 + (data[j][i,2:3] - data[1][i,2:3])**2 + (data[j][i,3:4] - data[1][i,4:5])**2)
                v = data[j][i,8:9]
                escapeV = np.sqrt((2 * sim.G * self.stars[0].m / d_1) + (2 * sim.G * self.stars[1].m / d_2))

                rnorm = data[j][i,4:5]
                rdot = data[j][i,5:8]
                r = data[j][i,1:4]
                a = (2*rnorm**-1 - np.dot(rdot,rdot))**-1
                e = np.sqrt(1 - (np.linalg.norm(np.cross(r,rdot)))**2 / a)

#                 if e >= 1:
#                     data[j][i,9:10] = 1
#                 else:
#                     data[j][i,9:10] = 0
                if v > escapeV:
                    data[j][i,9:10] = 1
                else:
                    data[j][i,9:10] = 0
                data[j][i, 10:11] = escapeV
            i += 1

        cols = "t x y z d vx vy vz v e ev".split()
        df_particles = [pd.DataFrame(data=data[i],columns=cols) for i in range(sim.N)]
        return df_particles


    def createSim(self,df_particles):
        sim = rebound.Simulation()
        sim.units = ("yr","AU","Msun")
        sim.G = 39.476926421373
        particle_index = 0
        for df in df_particles:
            if particle_index < len(self.stars):
                sim.add(m = self.stars[particle_index].m,
                        x = df['x'].iloc[-1],
                        y = df['y'].iloc[-1],
                        z = df['z'].iloc[-1],
                        vx = df['vx'].iloc[-1],
                        vy = df['vy'].iloc[-1],
                        vz = df['vz'].iloc[-1])
            else:
                sim.add(m = self.ast.m,
                        x = df['x'].iloc[-1],
                        y = df['y'].iloc[-1],
                        z = df['z'].iloc[-1],
                        vx = df['vx'].iloc[-1],
                        vy = df['vy'].iloc[-1],
                        vz = df['vz'].iloc[-1])
        return sim


    def consecutiveIntegrate(self,sim,t,no_of_sims,iterations_per_sim):
        t_gained = 0
        df_particles = self.integrateSim(sim,t_gained,t_gained+t/no_of_sims,iterations_per_sim,dt="linear")
        t += t/no_of_sims
        while t_gained < t:
            sim = self.createSim(df_particles)
            df = self.integrateSim(t_gained,t_gained+t/no_of_sims,iterations_per_sim,sim,dt="linear")
            for index in range(len(df_particles)):
                pd.concat([df_particles[index],df[index]])
            t_gained += t/no_of_sims
        return df_particles


    def mkpath(self,df_particles):
        """
        Creates path for data on filesystem

        Parameters
        ----------
        df_particles: list
            List of pandas DataFrames, each containing the t,x,y,z,d,vx,vy,vz,d
            of each star and asteroid for every integration iteration.

        Returns
        ----------
        path: String
            Name of the path created relative to ./
        """
        dir_name = ""
        for star in self.stars:
            dir_name += "Star({}M_{}AU)_".format(star.m,star.a)

        dir_name += "Asteroids({}AU_{}V_{}N)".format(self.ast.a,self.ast.vi,self.n_ast)

        t_f = int(df_particles[0]['t'].iloc[-1])
        iterations = df_particles[0].shape[0]
        sub_dir = "{}t_{}iterations".format(t_f,iterations)
        path = "{}/{}".format(dir_name,sub_dir)

        try:
            os.makedirs(path)
            print(dir_name + " created.")
        except FileExistsError:
            print(dir_name + " already exists.")
        finally:
            return path


    def finalPositionsVelocities(self,df_particles,printData=False,saveData=False):
        """
        Returns a dataframe containing the final positions and velocities
        of all the particles.

        Parameters
        ----------
        df_particles: list
            List of pandas DataFrames, each containing the t,x,y,z,d,vx,vy,vz,d
            of each star and asteroid for every integration iteration.
        printData: boolean
            Determines whether to print the data to the console.
        saveData: boolean
            Determines whether to save the data to the hard disk.

        Returns
        ----------
        df: pandas DataFrame
            DataFrame containing the final distances and velocities
            of each star and asteroid.
        """
        global G
        d_f = [particle['d'].iloc[-1] for particle in df_particles]
        v_f = [particle['v'].iloc[-1] for particle in df_particles]
        x_f = [particle['x'].iloc[-1] for particle in df_particles]
        y_f = [particle['y'].iloc[-1] for particle in df_particles]

        star_cols = ["star{}".format(i+1) for i in range(len(self.stars))]
        ast_cols = ["asteroid{}".format(i+1) for i in range(self.n_ast)]
        cols = star_cols + ast_cols

        df = pd.DataFrame(data=np.array([d_f,v_f, x_f, y_f]),columns=cols,index="d v x y".split()).T

        if printData:
            print(df.head())
        if saveData:
            path = self.mkpath(df_particles)
            df.to_csv("{}/finalPositionsVelocities.csv".format(path))
        return df


    def writeData(self,df_particles):
        """
        Creates folder and writes all the data to csv's.

        Parameters
        ----------
        df_particles: list
            List of pandas DataFrames, each containing the t,x,y,z,d,vx,vy,vz,d
            of each star and asteroid for every integration iteration.
        """
        path = self.mkpath(df_particles)
        for i in tqdm(range(len(df_particles))):
            if i < len(self.stars):
                df_particles[i].to_csv("{}/Star{}.csv".format(path,i+1))
            else:
                df_particles[i].to_csv("{}/Asteroid{}.csv".format(path,i+1-len(self.stars)))


    def plotAll(self,df_particles,dim,fname="plot",proj='xy',lim=10,savePlot=True):
        """
        Returns 2D or 3D plot of the given dataframe.

        Parameters
        ----------
        df_particles: list
            List of pandas DataFrames, each containing the t,x,y,z,d,vx,vy,vz,d
            of each star and asteroid for every integration iteration.
        fname: String
            Name of figure.
        dim: int
            2D or 3D (enter 2 or 3)
        savePlot: boolean
            Determines whether the plot is to be saved on the hard disk
        proj: String
            Determines what projection the 2 plot will be

        Returns
        ----------
        fig: matplotlib figure object
            2D or 3D plot of orbit

        """
        fig = plt.figure()
        if (dim!=2 and dim!=3):
            dim = int(input("Dimension must be 2 or 3: "))
            return self.plot(df_particles,fname,dim,proj,savePlot)
        elif (dim == 2):
            ax = fig.add_subplot(111)
            for i in range(len(df_particles)):
                if i < len(self.stars):
                    ax.plot(df_particles[i][proj[0]],
                            df_particles[i][proj[1]],"k")
                else:
                    ax.plot(df_particles[i][proj[0]],
                            df_particles[i][proj[1]])
        else:
            ax = fig.add_subplot(111, projection='3d')
            for i in range(len(df_particles)):
                if i < len(self.stars):
                    ax.plot(df_particles[i]['x'],
                            df_particles[i]['y'],
                            df_particles[i]['z'],"k")
                else:
                    ax.plot(df_particles[i]['x'],
                            df_particles[i]['y'],
                            df_particles[i]['z'])
            ax.set_zlim(-lim,lim)
            ax.set_zlabel("z (AU)")
        ax.set_xlim(-lim,lim)
        ax.set_ylim(-lim,lim)
        ax.set_xlabel("{} (AU)".format(proj[0]))
        ax.set_ylabel("{} (AU)".format(proj[1]))
        fname = "{}_{}D".format(fname,dim)
        plt.title(fname)
        path = self.mkpath(df_particles)
        if savePlot:
            plt.savefig("{}/{}.pdf".format(path,fname))
            plt.close(fig)
        return fig


    def plot(self,df_particles,particles=[0],fname="plot",lim=10,savePlot=True):

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for particle in particles:
            ax.plot(df_particles[particle]['x'],
                    df_particles[particle]['y'],
                    df_particles[particle]['z'])
        ax.set_zlim(-lim,lim)
        ax.set_xlim(-lim,lim)
        ax.set_ylim(-lim,lim)
        fname = "{}_{}D".format(fname, "3")
        plt.title(fname)
        path = self.mkpath(df_particles)
        if savePlot:
            plt.savefig("{}/{}.pdf".format(path,fname))
            plt.close(fig)
        return fig
