from rebound import Particle, StarSystem
from flask import Flask, flash, jsonify, redirect, session, render_template, request, session
from flask_session import Session


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
        return redirect("/Sim")
    else:
        return render_template("AstroSim.html")

@app.route("/Sim")
def Sim():
    return render_template("Sim.html")
