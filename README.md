# Adaptive Cruise Control using PID controller 
Adaptive cruise control with a PID (Proportional, Integral, Derivative) controller manages vehicle to vehicle distance based on the velocity or acceleration of leading vehicle.
The PID controller allows the ACC to automatically control the speed, throttle, and braking without driver intervention in certain conditions. 
The plant model is simplified and parameterized as a vehicle with two wheels, which uses Longitudinal Wheel models from Simulink Vehicle Dynamics Blockset. 

As it's shown in the below scope capture , in one hand we maintained a safe distance between the two cars, the distance is between 26m and 23m. On the other hand we reached a very low gap between the velocity of the leading( followed) and ego(following) vehicles.


![logs](https://user-images.githubusercontent.com/79842338/122574633-c3226a00-d047-11eb-8751-f38447a7e3e1.png)

To calculate the different physical parameter and Dynamics of the Vehcile we used the belwo specification/relation:

<br />Gravitational acceleration unit g : 9.81 [m/s2]
<br />Vehicle total mass M : 750 [kg]
<br />Front static load distribution ratio df : 0.6 [-]
<br />Rear static load distribution ration dr : 1 − df [-]
<br />Center of gravity height H : 0.55 [m]
<br />Wheelbase L : 2.6 [m]
<br />Front braking force distribution ratio bf : 0.8 [-]
<br />Rear braking force distribution ratio br : 1 − bf [-]
<br />Number of braking piston for one wheel n : 2 [-]
<br />Kinetic friction coefficient μ : 0.2 [-]
<br />Braking piston area A : 1.96 * 10−3[m2]
<br />Brake pad mean radius R : 0.1778 [m]

<br />Dynamic front weight can be calculated as below.
<br />Fz = ( M * df − M * a * H/L ) * g  Where a [m/s2] is the longitudinal acceleration.

<br />for the relation between the braking torque ( BrkTrq [Nm]) and braking piston pressure ( BrkPrs [Pa]), we used the equation below .
<br />BrkTrq = n * μ * BrkPrs * A * R



