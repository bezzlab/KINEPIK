from flask import Flask, flash, render_template, url_for, redirect, request, jsonify, json, Response
from KINEPIK_app.api_v1.protein_routes import protein_bp
from KINEPIK_app.api_v1.kinase_routes import kinase_bp
from KINEPIK_app.api_v1.inhibitor_routes import inhibitor_bp
from KINEPIK_app.api_v1.sif_routes import sif_bp
from KINEPIK_app.frontend.routes import routes_bp

# creting the app variable
app = Flask(__name__)

# registering the blueprints for the api
app.register_blueprint(protein_bp)
app.register_blueprint(kinase_bp)
app.register_blueprint(inhibitor_bp)
app.register_blueprint(sif_bp)

# registering the blueprints for the frontend
app.register_blueprint(routes_bp)

@app.route("/", methods = ["GET"])
def index():
    info = f"""
    <html>
        <body>
            <p> Arriving here shortly: KINEPIK! 😀😀 </p>
        </body>
    </html>
    """

    # return changes the format of the string to html that will be displayed to the user
    return Response(info, mimetype = "text/html")

# starting the server with debugging
if __name__ == "__main__":
    app.run(debug=True)