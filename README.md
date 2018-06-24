# Unscented Kalman Filter Project Starter Code
Term2-Project2: Unscented Kalman Filters

This project implements the Unscented Kalman filters in c++ for RADAR and LIDAR measurements. LIDAR data constitutes of x & y positions and RADAR data constitutes of radius, angle and radial velocity measurements.

The main difference in Unscented Kalman & Extended Kalman filters is, the UKF uses non-linear motion models to estimate the state of the object. Specifically it uses the CTRV (Constant Turn Rate and Velocity) model. The process logic remain the same with predict & update steps, except that, instead of using linear approximations like EKF, UKF uses sigma points to calculate the mean state and covariance. As a result, UKF is much better in predicting non-linear paths, with much lower RMSE.

This project works with the Udacity's simulator, which reads input data from the dataset file and the calculated RMSE is displayed on the simulator. To run the project, first start the simulator and then run the generated UKF binary from the project.

![RMSE](https://raw.githubusercontent.com/nitheeshkl/CarND-Unscented-Kalman-Filter-Project/master/rmse.png)