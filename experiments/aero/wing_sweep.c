#include <stdio.h>
#include <math.h>
#include <mujoco/mujoco.h>

static mjModel* M;

// Set angle of attack (pitch about Y, leading edge up) and wind speed,
// hold the wing at rest, and return aero force via qfrc_passive.
static void measure(double aoa_deg, double wind_x, double* lift, double* drag) {
  mjData* d = mj_makeData(M);
  double a = aoa_deg * M_PI / 180.0;
  d->qpos[3] = cos(a/2);    // quaternion, rotation about world Y
  d->qpos[4] = 0;
  d->qpos[5] = sin(a/2);
  d->qpos[6] = 0;
  M->opt.wind[0] = wind_x;  // override wind speed at runtime
  mju_zero(d->qvel, M->nv);
  mj_forward(M, d);
  *drag = d->qfrc_passive[0];   // Fx, along wind
  *lift = d->qfrc_passive[2];   // Fz, perpendicular (lift)
  mj_deleteData(d);
}

int main(int argc, char** argv) {
  char err[1000] = "";
  M = mj_loadXML(argc > 1 ? argv[1] : "/tmp/wing.xml", NULL, err, sizeof(err));
  if (!M) { printf("load error: %s\n", err); return 1; }

  printf("Wing: chord 0.4 m, span 2.0 m, area ~0.63 m^2 (ellipse).\n");

  printf("\n== Angle-of-attack sweep at 15 m/s ==\n");
  printf("  AoA     lift(N)    drag(N)    L/D\n");
  for (double aoa = 0; aoa <= 30.0001; aoa += 2.5) {
    double L, D;
    measure(aoa, 15.0, &L, &D);
    printf("  %4.1f   %8.2f   %8.3f   %6.2f\n", aoa, L, D, D != 0 ? L/D : 0);
  }

  printf("\n== Wind-speed sweep at AoA = 8 deg (lift ~ v^2) ==\n");
  printf("  wind(m/s)   lift(N)   lift/v^2\n");
  for (double v = 5; v <= 30.0001; v += 5) {
    double L, D;
    measure(8.0, v, &L, &D);
    printf("  %6.1f    %8.2f    %7.4f\n", v, L, L/(v*v));
  }

  // how fast must this wing fly to lift its own 2 kg (weight 19.6 N) at 8 deg?
  double L8, D8;
  measure(8.0, 1.0, &L8, &D8);            // lift coefficient at 1 m/s
  double v_needed = sqrt(19.62 / L8);     // since lift scales with v^2
  printf("\nTo lift its own weight (2 kg = 19.62 N) at 8 deg AoA: ~%.1f m/s airspeed.\n",
         v_needed);

  mj_deleteModel(M);
  return 0;
}
