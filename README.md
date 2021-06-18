# Adaptive Cruise Control using PID controller 
Adaptive cruise control with a PID (Proportional, Integral, Derivative) controller manages vehicle to vehicle distance based on the velocity or acceleration of leading vehicle.
The PID controller allows the ACC to automatically control the speed, throttle, and braking without driver intervention in certain conditions. 
The plant model is simplified and parameterized as a vehicle with two wheels, which uses Longitudinal Wheel models from Simulink Vehicle Dynamics Blockset. 

As it's shown in the below scope capture , in one hand we maintained a safe distance between the two cars, the distance is between 26m and 23m. On the other hand we reached a very low gap between the velocity of the leading( followed) and ego(following) vehicles.


![logs](https://user-images.githubusercontent.com/79842338/122574633-c3226a00-d047-11eb-8751-f38447a7e3e1.png)

