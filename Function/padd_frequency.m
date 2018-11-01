function fp=padd_frequency(f,padd)
df=max(diff(f));
fp_left=f(1)-df*[padd:-1:1];
fp_right=f(end)+df*[1:padd];

 fp=[fp_left f fp_right];
 