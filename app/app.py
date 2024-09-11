from flask import Flask, render_template, request
import math
import os
import matplotlib
matplotlib.use('Agg')  
import matplotlib.pyplot as plt
import numpy as np

app = Flask(__name__)
@app.route('/')
def home():
    return render_template('home.html')

@app.route('/wavelength', methods=['GET', 'POST'])
def index():
    result = None
    iterations = None
    celerity = None
    plot_url = None
    wave_type = None

    if request.method == 'POST':
        try:
            g = 9.81  # Acceleration due to gravity
            t = float(request.form['t'])
            d = float(request.form['d'])

            # Initial guess for wavelength
            input_l = (9.81 * t**2) / (2 * math.pi)
            cal_l = 0
            iteration_count = 0
            List_L = []
            List_it = []

            while True:
                cal_l = ((g * t**2) / (2 * math.pi)) * (math.tanh((2 * math.pi * d) / input_l))
                error = abs(cal_l - input_l)
                iteration_count += 1
                List_L.append(cal_l)
                List_it.append(iteration_count)

                if error < 0.001:
                    break
                input_l = cal_l

            result = "{:.2f}".format(cal_l)
            iterations = iteration_count
            celerity = "{:.2f}".format(cal_l / t)

            # Plot the convergence of L over iterations
            plt.figure()
            plt.plot(List_it, List_L, label='Calculated L')
            plt.xlabel('Iterations')
            plt.ylabel('Calculated L')
            plt.title('Convergence of L over Iterations')

            # Scatter plot and annotation
            plt.scatter(iterations, cal_l, color='red', zorder=5)
            plt.annotate(f'({iterations}, {result})', 
                         xy=(iterations, cal_l), 
                         xytext=(-40, 10), 
                         textcoords='offset points',
                         fontsize=10,
                         color='black',)  

            plot_filename = os.path.join(app.root_path, 'static', 'plot.png')
            plt.savefig(plot_filename)
            plt.close()

            plot_url = 'static/plot.png'

            # Determine wave type
            dl_ratio = d / cal_l

            if dl_ratio > 0.5:
                wave_type = "Deep Wave"
            elif dl_ratio > 0.05:
                wave_type = "Intermediate Wave"
            else:
                wave_type = "Shallow Wave"

        except ValueError as e:
            print(f"ValueError: {e}")
            result = "Error in input values"
        except Exception as e:
            print(f"Exception: {e}")
            result = "An error occurred"

    return render_template('index.html', result=result, iterations=iterations, 
                           celerity=celerity, plot_url=plot_url, wave_type=wave_type)

@app.route('/depth', methods=['GET', 'POST'])
def depth():
    result = None
    if request.method == 'POST':
        g = 9.8  
        t = float(request.form['t'])
        L = float(request.form['L'])

        if(2*math.pi*L/(9.81 * t**2)>1 or 2*math.pi*L/(g * t**2)< -1):
            result = "Math Error"

        else:
            result = (L/(2*math.pi))*math.atanh((2*math.pi*L/(g*t**2)))
            result = "{:.2f}".format(result)

    return render_template('depth.html', result=result)
  

@app.route('/vis', methods=['GET', 'POST'])
def vis():
  
    v = w = L = p = del_x = del_z = C = Cg = E = None
    displacement_plot_z_path = displacement_plot_path = pressure_plot_path = None
    a_x = a_z = acc_x_plot_path = acc_z_plot_path = None

    if request.method == 'POST':
     
        t = float(request.form['t'])
        d = float(request.form['d'])
        z = float(request.form['z'])
        h = float(request.form['h'])
        Q = float(request.form['Q'])
        Q = np.radians(Q)
        density = float(request.form['density'])
        g = 9.81
        gamma = density * g

        # Iterative calculation for wave length
        input_l = (9.81 * t**2) / (2 * math.pi)
        cal_l = 0
        while True:
            cal_l = ((g * t**2) / (2 * math.pi)) * (math.tanh((2 * math.pi * d) / input_l))
            error = abs(cal_l - input_l)
            if error < 0.001:
                break
            input_l = cal_l
        L = input_l
        k = 2 * math.pi / L
        sigma = 2 * math.pi / t
        C = L/t
        n = (1/2)*(1+(2*k*d/(math.sinh(2*k*d))))
        Cg = n*C
        E = (1/8)*(density)*(h**2)
    
        v = ((math.pi * h) / t) * math.cosh(k * (d + z)) * math.cos(Q) / math.sinh(k * d)
        w = -1 * ((math.pi * h) / t) * math.sinh(k * (d + z)) * math.sin(Q) / math.sinh(k * d)

        a_x = -1 * ((2 * math.pi * math.pi * h) / (t * t)) * math.cosh(k * (d + z)) * math.cos(Q) / math.sinh(k * d)
        a_z = -1 * ((2 * math.pi * math.pi * h) / (t * t)) * math.sinh(k * (d + z)) * math.sin(Q) / math.sinh(k * d)

        D = h * math.cosh(k * (d + z)) / (2 * math.sinh(k * d))
        B = h * math.sinh(k * (d + z)) / (2 * math.sinh(k * d))

        del_x = D * math.cos(Q)
        del_z = B * math.sin(Q)

        # Displacement plots
        x = np.linspace(0, 2 * np.pi, 100)
        x_degrees = np.degrees(x)  # Convert the x values (radians) to degrees
        y = ((math.pi * h) / t) * math.cosh(k * (d + z)) * np.cos(x) / math.sinh(k * d)

        plt.figure()
        plt.plot(x_degrees, y)
        plt.title('Horizontal paricle velocity vs Phase Angle (u vs θ)')
        plt.xlabel('Phase Angle (θ in degrees)')  # Update the label to indicate degrees
        plt.ylabel('horizontal particle velocity (u)')
        plt.grid(True)
        displacement_plot_path = os.path.join(app.root_path, 'static', 'displacement_plot.png')
        plt.savefig(displacement_plot_path)
        plt.close()


        y = -1 * ((math.pi * h) / t) * math.sinh(k * (d + z)) * np.sin(x) / math.sinh(k * d)

        plt.figure()
        plt.plot(x_degrees, y)  # Use x_degrees for the x-axis
        plt.title('Vertical paricle velocity vs Phase Angle (w vs θ)')
        plt.xlabel('Phase Angle (θ in degrees)')  # Update the label
        plt.ylabel('Vertical partical velocity (w)')
        plt.grid(True)
        displacement_plot_z_path = os.path.join(app.root_path, 'static', 'displacement_plot_z.png')
        plt.savefig(displacement_plot_z_path)
        plt.close()


        # Acceleration plots
        
        a1 = -1 * ((2 * math.pi * math.pi * h) / (t * t)) * math.cosh(k * (d + z))  / math.sinh(k * d)
        a2 = -1 * ((2 * math.pi * math.pi * h) / (t * t)) * math.sinh(k * (d + z))  / math.sinh(k * d)

        y1 = a1 * np.cos(x)
        

        plt.figure()
        plt.plot(np.degrees(x), y1)  # Convert x to degrees
        plt.title('Acceleration vs Phase Angle (ax vs θ)')
        plt.xlabel('Phase Angle (θ in degrees)')
        plt.ylabel('Acceleration ax')
        plt.grid(True)
        acc_x_plot_path = os.path.join(app.root_path, 'static', 'acc_x_plot.png')
        plt.savefig(acc_x_plot_path)
        plt.close()

        # Acceleration in the z-direction
       
        y = a_z * np.sin(x)
        plt.figure()
        plt.plot(np.degrees(x), y)  # Convert x to degrees
        plt.title('Acceleration vs Phase Angle (az vs θ)')
        plt.xlabel('Phase Angle (θ in degrees)')
        plt.ylabel('Acceleration az')
        plt.grid(True)
        acc_z_plot_path = os.path.join(app.root_path, 'static', 'acc_z_plot.png')
        plt.savefig(acc_z_plot_path)
        plt.close()


        # Pressure plot
        p = gamma * h / 2 * np.cosh(k * (d + z)) * np.sin(Q) / np.cosh(k * d) - gamma * z
        Q_vals = np.linspace(0, 2 * np.pi, 100)
        p_vals = gamma * h / 2 * np.cosh(k * (d + z)) * np.sin(Q_vals) / np.cosh(k * d) - gamma * z

        plt.figure()
        plt.plot(np.degrees(Q_vals), p_vals)  # Convert Q_vals to degrees
        plt.title('Pressure vs Phase Angle (p(θ))')
        plt.xlabel('Phase Angle (θ in degrees)')  # Update the label
        plt.ylabel('Pressure p(θ)')
        plt.grid(True)
        pressure_plot_path = os.path.join(app.root_path, 'static', 'pressure_plot.png')
        plt.savefig(pressure_plot_path)
        plt.close()

        def plot_velocity_field(H, T, d, lambda_wave, t=0.0, num_x=500, num_z=100):
            """
            Function to calculate and plot the horizontal water particle velocity field
            for a given wave height, period, depth, wavelength, and time.

            Parameters:
            H (float): Wave height in meters
            T (float): Wave period in seconds
            d (float): Water depth in meters
            lambda_wave (float): Wavelength in meters
            t (float, optional): Time in seconds. Default is 0.0.
            num_x (int, optional): Number of points in the x-direction. Default is 500.
            num_z (int, optional): Number of points in the z-direction. Default is 100.
            """
            k = 2 * np.pi / lambda_wave  # wave number (1/m)
            w = 2 * np.pi / T           # angular frequency (rad/s)

            # Create a meshgrid for x and z
            x_values = np.linspace(0, lambda_wave, num_x)  # horizontal position (m)
            z_values = np.linspace(0, d, num_z)            # depth from 0 (surface) to d (seabed)
            X, Z = np.meshgrid(x_values, z_values)

            # Calculate horizontal water particle velocity u(z, x, t)
            u_z_x = (H / T) * (np.cosh(k * (d - Z)) / np.sinh(k * d)) * np.cos(k * X - w * t)

            # Plotting the velocity field as a function of x and z
            plt.figure(figsize=(12, 8))
            plt.contourf(X, Z, u_z_x, cmap='viridis')
            plt.colorbar(label='Horizontal Velocity u(z, x, t) (m/s)')
            plt.title(f'Horizontal Water Particle Velocity u(z, x) at t = {t:.1f} s')
            plt.xlabel('Position x (m)')
            plt.ylabel('Depth z (m)')
            plt.gca().invert_yaxis()  # Invert the y-axis to have the surface at the top
            plt.grid(True)

            velocity_profile_path = os.path.join(app.root_path, 'static', 'velocity_profile.png')
            plt.savefig(velocity_profile_path)
            plt.close()

        plot_velocity_field(h, t, d, L, 0.0)

    return render_template(
        'vis.html',
        L=L, 
        displacement_plot_z=displacement_plot_z_path, 
        displacement_plot=displacement_plot_path, 
        pressure_plot=pressure_plot_path, 
        v=v, 
        w=w, 
        p=p, 
        del_x=del_x, 
        del_z=del_z,
        a_x=a_x,
        a_z=a_z,
        acc_x_plot_path=acc_x_plot_path,
        acc_z_plot_path=acc_z_plot_path,
        C=C,
        Cg=Cg,
        E=E
    )

@app.route('/trans', methods=['GET', 'POST'])
def trans():
    result = None
    shoaling_plot_path = None
    height_plot_path = None
    refraction_plot_path = None

    if request.method == 'POST':
        try:
            t = float(request.form.get('t', 0))
            d = float(request.form.get('d', 0))
            H0 = float(request.form.get('H0', 0))  # Deep water wave height
            alpha0 = float(request.form.get('alpha0', 0))  # Angle

            g = 9.81  # Gravitational constant
            input_l = (g * t**2) / (2 * math.pi)
            cal_l = 0
            while True:
                cal_l = ((g * t**2) / (2 * math.pi)) * (math.tanh((2 * math.pi * d) / input_l))
                error = abs(cal_l - input_l)
                if error < 0.001:
                    break
                input_l = cal_l
            L = input_l
            L0 = g * (t**2) / (math.pi * 2)

            # Function to calculate shoaling coefficient
            def calculate_shoaling_cof(T, d):
                g = 9.81  # Gravitational constant
                input_l = (g * t**2) / (2 * math.pi)
                cal_l = 0
                while True:
                    cal_l = ((g * t**2) / (2 * math.pi)) * (math.tanh((2 * math.pi * d) / input_l))
                    error = abs(cal_l - input_l)
                    if error < 0.001:
                        break
                    input_l = cal_l
                L = input_l
                k = 2 * math.pi / L  # Wave number
                n = 0.5 * (1 + (2 * k * d) / np.sinh(2 * k * d))  # Group velocity factor
                C0 = g * T / (2 * math.pi)  # Deep water wave speed
                C1 = L / T
                # Shoaling Coefficient
                Ks = math.sqrt(C0 / (2 * n * C1))
                return Ks
            
            Ks = calculate_shoaling_cof(t, d)
        
            # Generate and save the Shoaling Coefficient plot
            depths = np.linspace(0.1, d, 100)
            shoaling_coefficients = [calculate_shoaling_cof(t, dep) for dep in depths]
            plt.plot(depths, shoaling_coefficients)
            plt.xlabel('Water Depth (m)')
            plt.ylabel('Shoaling Coefficient (Ks)')
            plt.title('Shoaling Coefficient vs. Water Depth')
            plt.scatter(d, Ks, color='red', zorder=5)
            rounded_d = round(d, 2)
            rounded_ks = round(Ks, 2)
            plt.annotate(f'({rounded_d:.2f}, {rounded_ks:.2f})', 
             xy=(d, Ks), 
             xytext=(d-36, Ks+5),  # Offset for annotation text
             textcoords='offset points',
             fontsize=10,
             color='black')
            plt.grid(True)
            shoaling_plot_path = os.path.join(app.root_path, 'static', 'shoaling_coefficient_plot.png')
            plt.savefig(shoaling_plot_path)
            plt.close()

            # Convert alpha0 to radians for the refraction calculation
            alpha0_radians = np.radians(alpha0)

            # Function to calculate refraction coefficient
            def calculate_refraction_cof(alpha0_r, T, d):
                g = 9.81  # Gravitational constant
                input_l = (g * t**2) / (2 * math.pi)
                cal_l = 0
                while True:
                    cal_l = ((g * t**2) / (2 * math.pi)) * (math.tanh((2 * math.pi * d) / input_l))
                    error = abs(cal_l - input_l)
                    if error < 0.001:
                        break
                    input_l = cal_l
                L = input_l
                C = L / T
                C0 = g * T / (2 * math.pi)
                alpha = np.arcsin(C / C0 * np.sin(alpha0_r))
                Kr = np.sqrt(np.cos(alpha0_r) / np.cos(alpha))
                return Kr
            
            Kr = calculate_refraction_cof(alpha0_radians, t, d)
            h = Kr * Ks * H0

            # Generate and save the Refraction Coefficient plot
            alpha0_degrees = np.linspace(0, 90, 100)  # Range from 0 to 90 degrees
            alpha0_radians = np.radians(alpha0_degrees)  # Convert to radians

            # Calculate refraction coefficients for each alpha0
            refraction_coefficients = [calculate_refraction_cof(alpha, t, d) for alpha in alpha0_radians]

            # Generate the plot
            plt.plot(alpha0_degrees, refraction_coefficients)
            plt.scatter(alpha0, Kr, color='red', zorder=5)
            rounded_alpha0 = round(alpha0, 2)
            rounded_kr = round(Kr, 2)
            plt.annotate(f'({rounded_alpha0:.2f}, {rounded_kr:.2f})', 
             xy=(alpha0, Kr), 
             xytext=(alpha0-36, Kr+5),  # Offset for annotation text
             textcoords='offset points',
             fontsize=10,
             color='black')
            plt.xlabel('Angle α0 (degrees)')
            plt.ylabel('Refraction Coefficient (Kr)')
            plt.title('Refraction Coefficient vs. Angle')
            plt.grid(True)

            # Save the plot
            refraction_plot_path = os.path.join(app.root_path, 'static', 'refraction_coefficient_plot.png')
            plt.savefig(refraction_plot_path)
            plt.close()

            # Generate the Height vs Water Depth plot
            depths = np.linspace(0.1, d, 100)

            # Calculate heights for each depth (one set of heights for the entire range of depths)
            a = alpha0*math.pi / 180
            heights = [calculate_refraction_cof(a, t, dep) * calculate_shoaling_cof(t, dep) * H0 for dep in depths]

            # Clear the figure before plotting to avoid multiple plots stacking up
            plt.clf()

            # Plot the result
            plt.plot(depths, heights)
            plt.xlabel('Water Depth (m)')
            plt.ylabel('Height of Water Wave')
            plt.title('Transformed Wave Height vs Water Depth')
            plt.scatter(d, h, color='red', zorder=5)
            rounded_d = round(d, 2)
            rounded_h = round(h, 2)
            plt.annotate(f'({rounded_d:.2f}, {rounded_h:.2f})', 
             xy=(d, h), 
             xytext=(d-36, h+5),  # Offset for annotation text
             textcoords='offset points',
             fontsize=10,
             color='black')
            plt.grid(True)
            height_plot_path = os.path.join(app.root_path, 'static', 'height_vs_depth.png')
            plt.savefig(height_plot_path)
            plt.close()

            result = {
                "shoaling_coefficient": Ks,
                "refraction_coefficient": Kr,
                "Height": h
            }

        except Exception as e:
            result = {
                "error": str(e)
            }

    return render_template('trans.html', 
                           result=result, 
                           shoaling_plot_path='/static/shoaling_coefficient_plot.png', 
                           height_plot_path='/static/height_vs_depth.png', 
                           refraction_plot_path='/static/refraction_coefficient_plot.png')

if __name__ == '__main__':
    app.run(debug=True, host='127.0.0.1', port=5003)