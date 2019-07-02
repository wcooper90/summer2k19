# from . import Particle, StarSystem
from flask import Flask, flash, jsonify, redirect, session, render_template, request, session
from flask_session import Session


app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/error")
def error():
    return render_template("error.html")

@app.route("/Future")
def future():
    return render_template("future.html")
