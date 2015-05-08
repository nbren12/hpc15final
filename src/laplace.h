#ifdef __cplusplus
extern "C" {
#endif

void evolve_heat_equation_2d(double *x, int n, double dx,
			     int nt, double dt, int output_interval);

#ifdef __cplusplus
}
#endif
