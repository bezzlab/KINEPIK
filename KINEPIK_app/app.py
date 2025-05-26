from flask import Flask, flash, render_template, url_for, redirect, request, jsonify, json
from KINEPIK_app.api_v1.protein_routes import protein_bp
from KINEPIK_app.api_v1.kinase_routes import kinase_bp

# creting the app variable
app = Flask(__name__)

# registering all the blueprints
app.register_blueprint(protein_bp)
app.register_blueprint(kinase_bp)

# starting the server with debugging
if __name__ == "__main__":
    app.run(debug=True)