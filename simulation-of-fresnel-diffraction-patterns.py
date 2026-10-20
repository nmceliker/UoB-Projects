import scipy.integrate
import math as m
import numpy as np
import cmath as cm
import scipy
import matplotlib.pyplot as plt
from matplotlib import colormaps


def Fresnel1dreal(xp, x, k, z):  # integrand for the real part (cos)
    return np.cos(k / (2 * z) * (x - xp) ** 2)


def Fresnel1dimag(xp, x, k, z):  # integrand for the imaginary part (sin)
    return np.sin(k / (2 * z) * (x - xp) ** 2)


def Fresnel1dComplex(E_0, x1p, x2p):  # function integration real and imaginary parts of the integral
    result_real, _ = scipy.integrate.quad(
        lambda xp: Fresnel1dreal(xp, x, k, z), x1p, x2p
    )
    result_imag, _ = scipy.integrate.quad(
        lambda xp: Fresnel1dimag(xp, x, k, z), x1p, x2p
    )
    return (k * E_0) / (2 * np.pi * z) * result_real, (k * E_0) / (
        2 * np.pi * z
    ) * 1j * result_imag


# functions for the limits of integration (square aperture)
def xp1func(yp):
    if Input_shape == "c":
        return -np.sqrt(radius * radius - yp * yp)
    else:
        return yp1


def xp2func(yp):
    if Input_shape == "s":
        return yp2
    elif Input_shape == "r":
        return 3 * yp2
    elif Input_shape == "c":
        return np.sqrt(radius * radius - yp * yp)
    else:
        return yp2


def Fresnel2dreal(yp, xp, x, y, z):
    return np.cos(k / (2 * z) * ((x - xp) ** 2 + (y - yp) ** 2))


def Fresnel2dimag(yp, xp, x, y, z):
    return np.sin(k / (2 * z) * ((x - xp) ** 2 + (y - yp) ** 2))


def Fresnel2dComplex(x, y, z):
    realpart, real_error = scipy.integrate.dblquad(
        lambda xp, yp: Fresnel2dreal(xp, yp, x, y, z), yp1, yp2, xp1func, xp2func
    )
    imagpart, imag_error = scipy.integrate.dblquad(
        lambda xp, yp: Fresnel2dimag(xp, yp, x, y, z), yp1, yp2, xp1func, xp2func
    )
    return ((k * E_0) / (2 * np.pi * z)) * (realpart + 1j * imagpart)


MyInput = "0"
while MyInput != "q":
    MyInput = input('Enter a choice, "a", "b", or "q" to quit: ')
    print("You entered the choice: ", MyInput)

    if MyInput == "a":
        print("You have chosen part (a), 1-D Diffraction.")
        Input_awyn = input("Would you like to enter a specified aperture width (y/n): ")
        if Input_awyn == "y":
            Input_awyn = input("Enter a value for the aperture width (in metres): ")
            aperture_width = float(Input_awyn)
        elif Input_awyn == "n":
            aperture_width = 2e-5
            print(
                f"You have chosen not to specify an aperture width. Value previously set value of {aperture_width} m is gonna be used."
            )
        else:
            print("Invalid choice.")
        Input_z = input(
            "Would you like to enter a specified value for screen distance, z (y/n): "
        )
        if Input_z == "y":
            Input_z = input("Enter a value for z (in metres): ")
            z = float(Input_z)
        elif Input_z == "n":
            z = 0.02
            print(
                f"You have chosen not to specify the screen distance. Value previously set value of {z} m is gonna be used."
            )
        else:
            print("Invalid choice.")
        # variables and constants
        numpoints = 200
        xp1, xp2 = (
            -aperture_width / 2,
            aperture_width / 2,
        )  # upper and lower limits of int in metres (m)
        xp_values = np.linspace(xp1, xp2, numpoints)
        E_0 = 6.25
        lambd = 1e-6  # wavelength
        k = (2 * np.pi) / lambd
        x_values = np.linspace(-0.005, 0.005, numpoints)  # 0.005m

        # arrays for real and imaginary parts of E
        E_real = np.array([])
        E_imag = np.array([])

        for x in x_values:
            for xp in xp_values:
                real_part, imag_part = Fresnel1dComplex(E_0, xp1, xp2)

            E_real = np.append(E_real, real_part)
            E_imag = np.append(E_imag, imag_part)

        E_complex = (
            E_real + E_imag
        )  # array for E complex (with both real and imaginary parts)
        intensity_complex = (abs(E_complex)) ** 2  # intensity array

        plt.plot(x_values, intensity_complex)
        plt.title(f'1-D Diffraction Pattern with z = {z}')
        plt.xlabel("Screen coordinate (m)")
        plt.ylabel("Relative Intensity")
        plt.show()

    elif MyInput == "b":
        print("You have chosen part (b), 2-D Diffraction.")
        z = float(input('Enter a value for the screen distance, z (in metres):'))
        Input_shape = input(
            "Choose an aperture shape square, rectangle or circle (s/r/c):"
        )
        if Input_shape == "s":
            square_size_yn = input(
                "Would you like to specify the width of the square (y/n):"
            )
            if square_size_yn == "y":
                aperture_width = input(
                    "Enter a value for the area of the square (in metres):"
                )
                aperture_width = float(aperture_width)
            elif square_size_yn == "n":
                aperture_width = 2e-5

                print(
                    f"You have chosen not to specify the size, value of {aperture_width} will be used."
                )
            else:
                print("Invalid choice.")
            yp1 = -aperture_width / 2
            yp2 = aperture_width / 2
            # Variables
            numpoints = 100
            wavelength = 1e-6  # m
            k = 2 * np.pi / wavelength
            E_0 = 3

            # Grid of x, y values on the screen
            x_values = np.linspace(-5e-3, 5e-3, numpoints)  # in m
            y_values = np.linspace(-5e-3, 5e-3, numpoints)  # in m
            X, Y = np.meshgrid(x_values, y_values)

            intensity = np.zeros((numpoints, numpoints))
            for i, x in enumerate(x_values):
                for j, y in enumerate(y_values):
                    E = Fresnel2dComplex(x, y, z)
                    intensity[i, j] = abs(E) ** 2

            # Plot the intensity
            plt.figure()
            plt.imshow(intensity, extent=(-5, 5, -5, 5), origin="lower", cmap="hot")
            plt.colorbar(label="Intensity")
            plt.title(f"Diffraction Intensity for Square Aperture, z = {z} m")
            plt.xlabel("x (mm)")
            plt.ylabel("y (mm)")
            plt.show()

        elif Input_shape == "r":
            rectangle_size_yn = input(
                "Would you like to specify the width of the rectangle (y/n):"
            )
            if rectangle_size_yn == "y":
                aperture_width = input(
                    "Enter a value for the widdth of the rectangle (in metres):"
                )
                aperture_width = float(aperture_width)
            elif rectangle_size_yn == "n":
                aperture_width = 2e-5

                print(
                    f"You have chosen not to specify the size, value of {aperture_width} will be used."
                )
            else:
                print("Invalid choice.")
            yp1 = -aperture_width / 2
            yp2 = aperture_width / 2
            # Variables
            numpoints = 100
            wavelength = 1e-6  # m
            k = 2 * np.pi / wavelength
            E_0 = 3

            # Grid of x, y values on the screen
            x_values = np.linspace(-5e-3, 5e-3, numpoints)  # in m
            y_values = np.linspace(-5e-3, 5e-3, numpoints)  # in m
            X, Y = np.meshgrid(x_values, y_values)

            intensity = np.zeros((numpoints, numpoints))
            for i, x in enumerate(x_values):
                for j, y in enumerate(y_values):
                    E = Fresnel2dComplex(x, y, z)
                    intensity[i, j] = abs(E) ** 2

            # Plot the intensity
            plt.figure()
            plt.imshow(intensity, extent=(-5, 5, -5, 5), origin="lower", cmap="hot")
            plt.colorbar(label="Intensity")
            plt.title(f"Diffraction Intensity for Rectangle Aperture, z = {z} m")
            plt.xlabel("x (mm)")
            plt.ylabel("y (mm)")
            plt.show()

        elif Input_shape == "c":
            circle_size_yn = input(
                "Would you like to specify the radius of the circle (y/n):"
            )
            if circle_size_yn == "y":
                radius = input("Enter a value for the radius (in metres):")
                radius = float(radius)
            elif circle_size_yn == "n":
                radius = 2e-5
                print(
                    f"You have chosen not to specify the radius, {radius} m will be used."
                )
            else:
                print("Invalid choice.")
            yp1 = -radius
            yp2 = radius
            # Variables
            numpoints = 100
            wavelength = 1e-6  # m
            k = 2 * np.pi / wavelength
            E_0 = 3

            # Grid of x, y values on the screen
            x_values = np.linspace(-5e-3, 5e-3, numpoints)  # in m
            y_values = np.linspace(-5e-3, 5e-3, numpoints)  # in m
            X, Y = np.meshgrid(x_values, y_values)

            intensity = np.zeros((numpoints, numpoints))
            for i, x in enumerate(x_values):
                for j, y in enumerate(y_values):
                    E = Fresnel2dComplex(x, y, z)
                    intensity[i, j] = abs(E) ** 2

            # Plot the intensity
            plt.figure()
            plt.imshow(intensity, origin="lower", cmap="hot")
            plt.colorbar(label="Intensity")
            plt.title(f"Diffraction Intensity for Circle Aperture, z = {z} m")
            plt.xlabel("x (mm)")
            plt.ylabel("y (mm)")
            plt.show()

    elif MyInput != "q":
        print("This is not a valid choice")

print("You have chosen to finish - goodbye.")
