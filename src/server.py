from flask import Flask, Response, render_template, request, jsonify
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from simulation import Simulation
import numpy as np
import io
import webbrowser
import json

app = Flask(__name__)

current_simulation = None


@app.route("/", methods=["GET", "POST"])
def index():
    global current_simulation
    if request.method == "POST":
        data = json.loads(request.data)
        print(data)
        current_simulation = Simulation(
            bottle_volume=data["bottle_volume"],
            inside_pressure=data["bottle_pressure"],
            nozzle_area=data["nozzle_area"],
            launch_angle=data["launch_angle"],
            rocket_mass=data["rocket_mass"],
            fill_factor=data["fill_factor"],
            rocket_face_area=data["face_area"],
        )
        current_simulation.run()
        return jsonify(success=True)
    return render_template("index.html")


@app.route("/trajectory.png")
def trajectory():
    if current_simulation is None:
        return jsonify(success=False)
    fig = current_simulation.calc_trajectory(plot=True)
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype="image/png")


@app.route("/drag.png")
def drag():
    if current_simulation is None:
        return jsonify(success=False)
    fig = current_simulation.calc_rocket_drag(plot=True)
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype="image/png")


@app.route("/thrust.png")
def thrust():
    if current_simulation is None:
        return jsonify(success=False)
    fig = current_simulation.calc_rocket_thrust(plot=True)
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype="image/png")


@app.route("/water-speed.png")
def water_speed():
    if current_simulation is None:
        return jsonify(success=False)
    fig = current_simulation.calc_water_speed(plot=True)
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype="image/png")


@app.route("/water-vol.png")
def water_vol():
    if current_simulation is None:
        return jsonify(success=False)
    fig = current_simulation.calc_water_volume(plot=True)
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype="image/png")


@app.route("/animation")
def trajectory_animation():
    if current_simulation is None:
        return jsonify(success=False)
    times, points = current_simulation.calc_trajectory()
    return jsonify(times=times.tolist(), points=points.tolist())


if __name__ == "__main__":
    webbrowser.open("127.0.0.1:5000")
    app.run()
