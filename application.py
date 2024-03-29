from .rebound import Particle, StarSystem
from flask import Flask, flash, jsonify, redirect, session, render_template, request, session
from flask_session import Session
import numpy as np


app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/error")
def error():
    return render_template("error.html")

@app.route("/about")
def about():
    return render_template("about.html")

@app.route("/Future")
def future():
    return render_template("future.html")

@app.route("/AstroSim", methods=["GET", "POST"])
def AstroSim():
    if request.method == "POST":
        G = 39.476926421373
        ss = float(request.form.get("ss"))
        ns = float(request.form.get("ns"))
        na = int(request.form.get("na"))
        sa = float(request.form.get("sa")) * 0.000001
        time = float(request.form.get("time"))
        stars = [Particle(m=1,x=9),Particle(m=1,x=-9)]
        stars[0].vy = np.sqrt(G*stars[1].m/(2*(stars[0].a+stars[1].a)))
        stars[1].vy = -np.sqrt(G*stars[0].m/(2*(stars[0].a+stars[1].a)))
        ast = Particle(m=sa,x=10)
        sys = StarSystem(stars,ast,n_ast=na,seed=1)
        t_start = 0
        t_stop = time
        iterations = 100
        sim = sys.initiateSim()
        df_particles = sys.integrateSim(t_start,t_stop,iterations,sim,dt="linear")
        sys.writeData(df_particles)
        df = sys.finalPositionsVelocities(df_particles,printData=True,saveData=True)
        sys.plotAll(df_particles,dim=3,lim=25)
        # maybe move this function into the rebound file. Figure out how to display the pdf without having to
        # save it somewhere, as well as links to the downloadable csv files for particles.
        info = sys.mkpath(df_particles)
        return render_template("Sim.html", info=info)
    else:
        return render_template("AstroSim.html")
