function plotInitialGuess(col_pts, x_guess, h_guess, v_guess, gamma_guess, m_guess, U_LG_guess)
% this function takes initial guesses for states and controls
% and plot them in several figures

figure
plot(x_guess/1000, h_guess/1000);
xlabel("x [km]")
ylabel("h [km]")

figure
plot(col_pts, x_guess / 1000);
xlabel("time [sec]")
ylabel("x [km]")

figure
plot(col_pts, h_guess / 1000);
xlabel("time [sec]")
ylabel("h [km]")

figure
plot(col_pts, v_guess);
xlabel("time [sec]")
ylabel("v [m/s]")

figure
plot(col_pts, gamma_guess * 180/pi);
xlabel("time [sec]")
ylabel("\gamma [deg]")

figure
plot(col_pts, m_guess);
xlabel("time [sec]")
ylabel("mass [kg]")

figure
plot(col_pts, U_LG_guess * 180/pi);
xlabel("time [sec]")
ylabel("\alpha [deg]")


end