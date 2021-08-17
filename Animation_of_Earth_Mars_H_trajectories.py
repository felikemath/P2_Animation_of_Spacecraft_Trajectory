import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.patches import Arc
from matplotlib import image
import matplotlib as mpl
import math


def generate_animation_points(R_Earth=1.0, Angle_Earth_0_deg = 0, R_Mars=1.5237, Angle_Mars_0_deg = None, num_points=2070):
    # Hohmann transfer orbit: A is half length of the semi-major axis, B is half length of the semi-minor axis,
    # C is half of the distance between the two foci
    A = (R_Mars + R_Earth) / 2
    C = (R_Mars - R_Earth) / 2
    B = math.sqrt(A * A - C * C)

    # Kepler's Third Law: the squares of the orbital periods of the planets are directly proportional to
    # the cubes of the semi major axes of their orbits,
    # Orbit period ratio of Earth : Hohmann : Mars = sqrt(R_Earth^3 : A^3 : R_Mars^3)
    # Therefore, During the time length (T) of half period of the Hohmann Orbit (for pi), the angles (in radians)
    # that Earth and Mars travel are, respectively
    angle_range_Earth = math.pi * math.sqrt(A**3 / R_Earth**3)
    angle_range_Mars = math.pi * math.sqrt(A ** 3 / R_Mars ** 3)

    Angle_Earth_0 = Angle_Earth_0_deg * math.pi / 180   # assume Earth starts from 0-degree (3-O'Clock direction)

    if Angle_Mars_0_deg is None:  # if not specified, calculate it
        Angle_Mars_0_deg = 180 * (1-math.sqrt(A**3 / R_Mars**3))
    Angle_Mars_0 = Angle_Mars_0_deg * math.pi / 180

    angle_start_Earth, angle_end_Earth = Angle_Earth_0, Angle_Earth_0+angle_range_Earth
    angle_start_Mars, angle_end_Mars = Angle_Mars_0, Angle_Mars_0+angle_range_Mars

    # Assume both Earth and Mars are orbiting around the sun in approximately circular orbits with constant angular speeds.
    # Create num_points sampling points (evenly sampled in TIME-domain over the duration of T)
    angles_Earth = np.linspace(angle_start_Earth, angle_end_Earth, num_points)
    Xloc_Earth = R_Earth * np.cos(angles_Earth)
    Yloc_Earth = R_Earth * np.sin(angles_Earth)

    angles_Mars = np.linspace(angle_start_Mars, angle_end_Mars, num_points)
    Xloc_Mars = R_Mars * np.cos(angles_Mars)
    Yloc_Mars = R_Mars * np.sin(angles_Mars)

    # Please note that, for the Hohmann transfer orbit, the speed is keep changing.
    # It's more convenient to calculate ellipse in polar coordinates
    e = C / A  # e is the eccentricity of the elliptical orbit
    theta = 0.0
    Xloc_Hohmann = [R_Earth]
    Yloc_Hohmann = [0.0]
    TravelLen_Hohmann = [0.0]

    for i in range(1, num_points):
        numerator = (math.pi * B / num_points / A) * (1 + e * math.cos(theta))**2
        denominator = (1 - e**2)**2 + (math.pi * B / num_points / A) * (1 + e * math.cos(theta)) * e * math.sin(theta)
        theta_step = numerator / denominator  # in unit of radian
        old_r = A * (1 - e ** 2) / (1 + e * math.cos(theta))
        theta += theta_step  # in unit of radian
        new_r = A * (1 - e ** 2) / (1 + e * math.cos(theta))
        newX = new_r * math.cos(theta)
        newY = new_r * math.sin(theta)

        Xloc_Hohmann.append(newX)
        Yloc_Hohmann.append(newY)

        # approximate travel length
        travelLen = TravelLen_Hohmann[i-1] + theta_step * (old_r+new_r) / 2  # in unit of AU
        TravelLen_Hohmann.append(travelLen)

        todalDays = int(365 * angle_range_Earth / (2*math.pi))  # approximately

    return (Xloc_Earth, Yloc_Earth, Xloc_Mars, Yloc_Mars, Xloc_Hohmann, Yloc_Hohmann, TravelLen_Hohmann, todalDays)


# initialization function for the animation: reset the x- and y-coordinate lists and the texts
def init():
    line_Hohmann.set_data([], [])
    line_Earth.set_data([], [])
    line_Mars.set_data([], [])

    text_day.set_text('')
    text_dist_2_Earth.set_text('')
    text_dist_2_Mars.set_text('')
    text_dist_2_Sun.set_text('')
    text_spacecraft_speed.set_text('')
    text_spacecraft_travellen.set_text('')

    return line_Earth, line_Mars, line_Hohmann, imgEarthPlot, imgMarsPlot, imgSatPlot, \
           text_day, text_dist_2_Mars, text_dist_2_Sun, text_dist_2_Earth, text_spacecraft_speed, text_spacecraft_travellen


# animation step function
def animate(i):
    if i == 0:  # clean up the previous plots
        xlist_Hohmann[:] = []
        ylist_Hohmann[:] = []
        xlist_Earth[:] = []
        ylist_Earth[:] = []
        xlist_Mars[:] = []
        ylist_Mars[:] = []

    x_Hohmann, y_Hohmann = Xloc_Hohmann[i], Yloc_Hohmann[i]
    xlist_Hohmann.append(x_Hohmann)
    ylist_Hohmann.append(y_Hohmann)
    line_Hohmann.set_data(xlist_Hohmann, ylist_Hohmann)

    x_Earth, y_Earth = Xloc_Earth[i], Yloc_Earth[i]
    xlist_Earth.append(x_Earth)
    ylist_Earth.append(y_Earth)
    line_Earth.set_data(xlist_Earth, ylist_Earth)


    x_Mars, y_Mars = Xloc_Mars[i], Yloc_Mars[i]
    xlist_Mars.append(x_Mars)
    ylist_Mars.append(y_Mars)
    line_Mars.set_data(xlist_Mars, ylist_Mars)

    imgMarsPlot.set_transform(mpl.transforms.Affine2D().translate(x_Mars - imgInitLocX, y_Mars - imgInitLocY) + ax.transData)
    imgEarthPlot.set_transform(mpl.transforms.Affine2D().translate(x_Earth-imgInitLocX, y_Earth-imgInitLocY) + ax.transData)
    imgSatPlot.set_transform(mpl.transforms.Affine2D().translate(x_Hohmann - imgInitLocX, y_Hohmann - imgInitLocY) + ax.transData)

    day = int(i / num_points * totalTravelDays + 0.5)
    text_day.set_text('Earth Days = {:3d}'.format(day))
    travel_length_Hohmann = int(TravelLen_Hohmann[i] * AU)

    dist_2_Earth = int(math.sqrt((x_Hohmann-x_Earth)**2+(y_Hohmann-y_Earth)**2) * AU)
    text_dist_2_Earth.set_text('Distance to Earth = {:,d} km'.format(dist_2_Earth))
    dist_2_Mars = int(math.sqrt((x_Hohmann - x_Mars) ** 2 + (y_Hohmann - y_Mars) ** 2) * AU)
    text_dist_2_Mars.set_text('Distance to Mars = {:,d} km'.format(dist_2_Mars))
    dist_2_Sun = int(math.sqrt((x_Hohmann) ** 2 + (y_Hohmann) ** 2) * AU)
    text_dist_2_Sun.set_text('Distance to Sun = {:,d} km'.format(dist_2_Sun))
    uSun = 1.32712440018E20 # m^3/s^2
    A = (R_Mars + R_Earth) / 2 * AU
    speed = math.sqrt(uSun * (2/dist_2_Sun/1000 - 1/A/1000) )  # m/s
    speed /= 1000  # km/s
    text_spacecraft_speed.set_text('Spacecraft\'s speed = {:.3f} km/s'.format(speed))
    text_spacecraft_travellen.set_text('Spacecraft has traveled $ \\approx $ {:,d} km'.format(travel_length_Hohmann))

    return line_Earth, line_Mars, line_Hohmann, imgEarthPlot, imgMarsPlot, imgSatPlot, \
           text_day, text_dist_2_Mars, text_dist_2_Sun, text_dist_2_Earth, text_spacecraft_speed, text_spacecraft_travellen

# set parameters
AU = 149598000.0  # the constant of Astronomical Unit, in km
R_Earth = 1.0  # Radius of Earth in AU
Angle_Earth_0_deg = 0
R_Mars=1.5237  # Radius of Mars in AU
Angle_Mars_0_deg = 44.287  # Mars is ahead of the Earth by this angle in degree
num_points = 2070  # number of sampling points for the animation

# 1. generate the animation data for the motions of Earth, Mars and Spacecraft on Hohmann Transfer Orbit
Xloc_Earth, Yloc_Earth, Xloc_Mars, Yloc_Mars, Xloc_Hohmann, Yloc_Hohmann, TravelLen_Hohmann, totalTravelDays = \
    generate_animation_points(R_Earth=R_Earth,
                              Angle_Earth_0_deg=Angle_Earth_0_deg,
                              R_Mars=R_Mars,
                              Angle_Mars_0_deg=Angle_Mars_0_deg,
                              num_points=num_points)

# 2. generate the static data of the full orbits of Earth and Mars ( to be plotted as dashed lines)
angles_0_to_2pi = np.linspace(0, 2 * math.pi, 200)
X_full_Earth = R_Earth * np.cos(angles_0_to_2pi)
Y_full_Earth = R_Earth * np.sin(angles_0_to_2pi)
X_full_Mars = R_Mars * np.cos(angles_0_to_2pi)
Y_full_Mars = R_Mars * np.sin(angles_0_to_2pi)


# 3. Set up the figure and ax, plot the basic elements
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
plt.subplots_adjust(left=0.08, bottom=0.08, right=0.92, top=0.92)
plt.title("Animation of Earth, Mars and Spacecraft Motions", fontsize=24)
plt.axhline(y=0.0, color='gray', linestyle='--', linewidth=0.5)
plt.axvline(x=0.0, color='gray', linestyle='--', linewidth=0.5)
plt.xlim(-2, 2)
plt.ylim(-2, 2)
ax.set_aspect('equal', adjustable='box')
plt.plot(X_full_Earth, Y_full_Earth, color="grey", linewidth=0.5, linestyle="--")
plt.plot(X_full_Mars, Y_full_Mars, color="grey", linewidth=0.5, linestyle="--")
plt.plot(Xloc_Hohmann, Yloc_Hohmann, color="grey", linewidth=0.5, linestyle="--")
plt.text(0.06, -0.12, 'Sun', color='white')
plt.annotate('Earth orbit', xy=(-0.9, -0.43), xycoords='data',
             xytext=(-0.75, -0.25), textcoords='data', fontsize=10, color='white',
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2", color='grey'))
plt.annotate('Mars orbit', xy=(-0.9, -1.23), xycoords='data',
             xytext=(-0.7, -1.2), textcoords='data', fontsize=10, color='white',
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2", color='grey'))
plt.annotate('Hohmann\ntransfer\norbit', xy=(0.75, 0.73), xycoords='data',
             xytext=(0.9, 0.8), textcoords='data', fontsize=10, color='white', va='top',
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2", color='grey'))
plt.annotate('Hohmann \nperihelion\n(Launch)', xy=(0.99, -0.01), xycoords='data',
             xytext=(0.45, -0.20), textcoords='data', fontsize=10, color='white', va='top',
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2", color='grey'))
plt.annotate('Hohmann \naphelion\n(Landing)', xy=(-1.52, 0), xycoords='data',
             xytext=(-1.92, 0.35), textcoords='data', fontsize=10, color='white', va='top',
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2", color='grey'))

ax.set_xticks([])
ax.set_yticks([])

# x- and y-coordinate lists of the spacecraft, Earth and Mars to be used during animation
xlist_Hohmann = []
ylist_Hohmann = []
xlist_Earth = []
ylist_Earth = []
xlist_Mars = []
ylist_Mars = []

# generate the plotting artists for the trajectories
# point_Hohmann, = ax.plot([], [], 'go', lw=0.5, markersize=4)
line_Hohmann, = ax.plot([], [], '-g', lw=1)
# point_Earth, = ax.plot([], [], 'bo', lw=0.5, markersize=4)
line_Earth, = ax.plot([], [], '-', lw=1, color='#3280FF')
# point_Mars, = ax.plot([], [], 'ro', lw=0.5, markersize=4)
line_Mars, = ax.plot([], [], '-r', lw=1)

# generate the plotting artists for the information texts
y_display_loc = 1.0
text_day = ax.text(0.02, y_display_loc - 0.03, '', transform=ax.transAxes, color='white')
text_spacecraft_speed = ax.text(0.02, y_display_loc - 2*0.03, '', transform=ax.transAxes, color="lightgreen")
text_spacecraft_travellen = ax.text(0.02, y_display_loc - 3*0.03, '', transform=ax.transAxes, color="lightgreen")
text_dist_2_Earth = ax.text(0.6, y_display_loc - 0.03, '', transform=ax.transAxes, color='white')
text_dist_2_Mars = ax.text(0.6, y_display_loc - 2*0.03, '', transform=ax.transAxes, color='white')
text_dist_2_Sun = ax.text(0.6, y_display_loc - 3*0.03, '', transform=ax.transAxes, color='white')

imgUniverse = image.imread('universe.png')
imgUniversePlot = ax.imshow(imgUniverse, extent= (-2, 2, -2, 2))
imgInitLocX = 10.0
imgInitLocY = 0.0
ext = (imgInitLocX-0.05, imgInitLocX+0.05, imgInitLocY-0.05, imgInitLocY+0.05)
imgEarth = image.imread('Earth32x32.png')
imgEarthPlot = ax.imshow(imgEarth, extent= ext)
imgMars  = image.imread('Mars32x32.png')
imgMarsPlot = ax.imshow(imgMars, extent= ext)
imgSat  = image.imread('satellite32x32.png')
imgSatPlot = ax.imshow(imgSat, extent= ext)
imgSun = image.imread('Sun.png')
imgSunPlot = ax.imshow(imgSun, extent= (-0.093, 0.107, -0.1, 0.1))




def main():
    # call the animator.
    intervalms = 1  # this means 1 ms per frame for display. This has no impact on anim.save(), which is controlled by fps.
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=num_points, interval=intervalms, blit=True)
    # save the animation as an mp4.  This requires ffmpeg to be installed.
    anim.save(r'output/Animation_Of_Planet_Orbits.mp4', fps=30, bitrate=1800, extra_args=['-vcodec', 'h264'])
    plt.show()


if __name__=='__main__':
    main()