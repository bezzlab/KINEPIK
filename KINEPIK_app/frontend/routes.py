from flask import Flask, flash, render_template, url_for, redirect, request, jsonify, json,Blueprint

routes_bp = Blueprint("routes_bp", __name__, url_prefix="/testing")

@routes_bp.route("/", methods = ["GET"])
def index():
    return render_template("index.html")

@routes_bp.route("/api", methods = ["GET"])
def api():
    return render_template("api.html")

@routes_bp.route("/about", methods = ["GET"])
def about():
    return render_template("about.html")
