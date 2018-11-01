function yp=padd_signal(y,padd)
 yp=[y(1)*ones(1,padd) y y(end)*ones(1,padd)];
 