import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate


# PARAMETERS
R = 8.314
INITIAL_TEMP = 293
ATM = 1e5
LIQUID_DENSITY = 1000  # Kg / m^3
AIR_DENSITY = 1.225  # Kg / m^3
G = -9.8
G_VEC = np.array([0, G])
GAMMA = 1.4
AIR_DRAG_COEFF = 0.75


class Simulation:
    def __init__(
        self,
        *,
        bottle_volume=None,
        inside_pressure=None,
        nozzle_area=None,
        launch_angle=None,
        rocket_mass=None,
        fill_factor=None,
        rocket_face_area=None,
    ):

        self.Vb = bottle_volume
        self.P0 = inside_pressure
        self.NA = nozzle_area
        self.launch_angle = launch_angle
        self.M = rocket_mass
        self.V0 = bottle_volume * (1 - fill_factor)
        self.V_liquid = bottle_volume * (fill_factor)
        self.FA = rocket_face_area
        self.K = inside_pressure * (bottle_volume * (1 - fill_factor)) ** GAMMA
        self.initial_speed = np.sqrt(2 * (inside_pressure - ATM) / LIQUID_DENSITY)
        self.drag_coeff = AIR_DRAG_COEFF * 0.5 * AIR_DENSITY * self.FA

        if self.V_liquid > 1:
            raise ValueError("Fill factor too high")

        self._run = False

    def run(self):
        if self._run:
            raise Exception("Already run simulation")

        self._calc_launch()
        self._calc_flight()

        self._run = True

    def _calc_launch(self):
        def water_exit_acc(t, vel):
            A = -(GAMMA * self.K * self.NA) / LIQUID_DENSITY
            B = 0.5 * LIQUID_DENSITY * vel ** 2 + ATM
            return A * (B / self.K) ** ((GAMMA + 1) / GAMMA)

        def find_limit(vel_func):
            expected_area = self.V_liquid / self.NA
            for limit in np.linspace(0, 1, 1000):
                if integrate.quad(vel_func, 0, limit)[0] > expected_area:
                    return limit
            return np.inf

        r_water_vel = integrate.solve_ivp(
            water_exit_acc, (0, 1), np.array([self.initial_speed]), dense_output=True
        )

        self.water_vel = np.vectorize(r_water_vel.sol)
        self.launch_time = find_limit(self.water_vel)

        total_mass = self.M + self.V_liquid * LIQUID_DENSITY

        def launch_acceleration(t, v):
            vm = np.linalg.norm(v)
            um = self.water_vel(t)
            mass = (
                total_mass
                - integrate.quad(self.water_vel, 0, t)[0] * self.NA * LIQUID_DENSITY
            )
            rocket_acc = (
                um ** 2 * self.NA * LIQUID_DENSITY - (vm ** 2 * self.drag_coeff)
            ) / mass
            if vm != 0:
                rocket_dir = v / vm
            else:
                rocket_dir = [
                    np.cos(self.launch_angle / 180 * np.pi),
                    np.sin(self.launch_angle / 180 * np.pi),
                ]
            return rocket_acc * rocket_dir + G_VEC

        self.instant_impulse = (
            lambda t: self.water_vel(t) ** 2 * self.NA * LIQUID_DENSITY
        )
        self.launch_impulse = integrate.quad(self.instant_impulse, 0, self.launch_time)[
            0
        ]

        r_launch_vel = integrate.solve_ivp(
            launch_acceleration,
            (0, self.launch_time),
            np.array([0, 0]),
            dense_output=True,
        )

        self.launch_vel = r_launch_vel.sol
        self.launch_vel_vec = np.vectorize(lambda t, axis: r_launch_vel.sol(t)[axis])
        self.launch_vel_mod = np.vectorize(
            lambda t: np.linalg.norm(r_launch_vel.sol(t))
        )

        t = np.linspace(0, self.launch_time, 25)

        dx = self.launch_vel_vec(t, 0)
        dy = self.launch_vel_vec(t, 1)

        x = integrate.cumtrapz(dx, t, initial=0)[-1]
        y = integrate.cumtrapz(dy, t, initial=0)[-1]

        self.post_launch_pos = np.array([x, y])
        if self.post_launch_pos[1] <= 0:
            print("[Warning]: Values lead to a rocket that never leaves the floor")

    def _calc_flight(self):
        def drag_acceleration(t, vel):
            vm = np.linalg.norm(vel)
            return -(vm * self.drag_coeff / self.M) * vel + G_VEC

        initial_velocity = self.launch_vel(self.launch_time)

        r_vel_traject = integrate.solve_ivp(
            drag_acceleration, (0, 100), initial_velocity, dense_output=True
        )
        self.flight_vel_vec = np.vectorize(lambda t, axis: r_vel_traject.sol(t)[axis])
        self.flight_vel_mod = np.vectorize(
            lambda t: np.linalg.norm(r_vel_traject.sol(t))
        )

        t = np.linspace(0, 10, 1000)

        dy = self.flight_vel_vec(t, 1)
        y = integrate.cumtrapz(dy, t, initial=0) + self.post_launch_pos[1]

        self.flight_time = t[np.argmax(y < 0)]

    def calc_water_speed(self, plot=False):
        if not self._run:
            raise Exception("no values calculated yet")

        t = np.linspace(0, self.launch_time, 50)
        water_speed = self.water_vel(t)

        if plot:
            fig = plt.figure(figsize=(5, 5))
            plot = fig.add_subplot(111)
            plot.plot(t, water_speed)
            plot.set_xlabel("Time since launch (s)")
            plot.set_ylabel("Water Speed at nozzle (m/s)")
            plot.set_title("Water exit velocity over time")
            plot.grid(linestyle="--")
            return fig
        else:
            return t, water_speed

    def calc_water_volume(self, plot=False):
        if not self._run:
            raise Exception("no values calculated yet")

        t = np.linspace(0, self.launch_time, 1000)
        volume = (
            self.V_liquid
            - np.cumsum(
                np.hstack(
                    [[0], self.water_vel(t[:999]) * self.NA * self.launch_time / 1000]
                )
            )
        ) * 1000

        if plot:
            fig = plt.figure(figsize=(5, 5))
            plot = fig.add_subplot(111)
            plot.plot(t, volume)
            plot.set_xlabel("Time since launch (s)")
            plot.set_ylabel("Water Volume (litres)")
            plot.set_title("Water volume over time")
            plot.grid(linestyle="--")
            return fig
        else:
            return t, volume

    def calc_rocket_thrust(self, plot=False):
        if not self._run:
            raise Exception("no values calculated yet")

        t = np.linspace(0, self.launch_time, 50)
        water_thrust = self.instant_impulse(t)

        if plot:
            fig = plt.figure(figsize=(5, 5))
            plot = fig.add_subplot(111)
            plot.plot(t, water_thrust)
            plot.set_xlabel("Time since launch (s)")
            plot.set_ylabel("Thrust (n)")
            plot.set_title("Thrust over time")
            plot.grid(linestyle="--")
            return fig
        else:
            return t, water_thrust

    def calc_rocket_drag(self, plot=False):
        if not self._run:
            raise Exception("no values calculated yet")

        t1 = np.linspace(0, self.launch_time, 50)
        t2 = np.linspace(0, self.flight_time, 50)
        drag1 = self.launch_vel_mod(t1) ** 2 * self.drag_coeff
        drag2 = self.flight_vel_mod(t2) ** 2 * self.drag_coeff

        t = np.hstack([t1, t2 + self.launch_time])
        drag = np.hstack([drag1, drag2])

        if plot:
            fig = plt.figure(figsize=(5, 5))
            plot = fig.add_subplot(111)
            plot.plot(t, drag)
            plot.axvline(x=self.launch_time, color="red", linestyle="--")
            plot.set_xlabel("Time since launch (s)")
            plot.set_ylabel("Drag (n)")
            plot.set_title("Drag over time")
            plot.grid(linestyle="--")
            return fig
        else:
            return t, drag

    def calc_trajectory(self, plot=False):

        if not self._run:
            raise Exception("no values calculated yet")

        t_launch = np.linspace(0, self.launch_time, 50)
        t_flight = np.linspace(0, self.flight_time, 200)

        dx_launch = self.launch_vel_vec(t_launch, 0)
        dy_launch = self.launch_vel_vec(t_launch, 1)
        dx_flight = self.flight_vel_vec(t_flight, 0)
        dy_flight = self.flight_vel_vec(t_flight, 1)

        x_launch = integrate.cumtrapz(dx_launch, t_launch, initial=0)
        y_launch = integrate.cumtrapz(dy_launch, t_launch, initial=0)

        x_flight = (
            integrate.cumtrapz(dx_flight, t_flight, initial=0) + self.post_launch_pos[0]
        )
        y_flight = (
            integrate.cumtrapz(dy_flight, t_flight, initial=0) + self.post_launch_pos[1]
        )

        t = np.hstack([t_launch, t_flight + self.launch_time])
        x = np.hstack([x_launch, x_flight])
        y = np.hstack([y_launch, y_flight])

        if plot:
            fig = plt.figure(figsize=(5, 5))
            plot = fig.add_subplot(111)
            plot.plot(x_launch, y_launch, color="red")
            plot.plot(x_flight, y_flight, color="C0")
            plot.set_xlabel("Distace (m)")
            plot.set_ylabel("Height (m)")
            plot.set_title("Trajectory")
            plot.grid(linestyle="--")
            plot.axis("equal")
            return fig
        else:
            return t, np.vstack([x, y]).T
