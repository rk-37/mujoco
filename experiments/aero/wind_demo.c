#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mujoco/mujoco.h>

// Load model, set the plate orientation (rotation about Y in radians),
// step once with the body at rest, and report the passive (fluid) force.
static void run(const char* xml, double pitch_rad, const char* label) {
  char err[1000] = "";
  mjModel* m = mj_loadXML(xml, NULL, err, sizeof(err));
  if (!m) { printf("load error: %s\n", err); return; }
  mjData* d = mj_makeData(m);

  // orient the plate: quaternion for rotation about world Y axis
  d->qpos[3] = cos(pitch_rad/2);   // w
  d->qpos[4] = 0;                  // x
  d->qpos[5] = sin(pitch_rad/2);   // y
  d->qpos[6] = 0;                  // z
  // body starts at rest; wind (5 m/s along +x) is the only relative flow
  mju_zero(d->qvel, m->nv);

  // compute forces without integrating motion away
  mj_forward(m, d);

  // free joint: qfrc_passive[0:3] = force (world frame), [3:6] = torque (body frame)
  double fx = d->qfrc_passive[0], fy = d->qfrc_passive[1], fz = d->qfrc_passive[2];
  double tx = d->qfrc_passive[3], ty = d->qfrc_passive[4], tz = d->qfrc_passive[5];
  double fmag = sqrt(fx*fx + fy*fy + fz*fz);

  printf("%-22s  drag(Fx)=%+.5f N  lift(Fz)=%+.5f N  |F|=%.5f N  torque_y=%+.5f N·m\n",
         label, fx, fz, fmag, ty);

  mj_deleteData(d);
  mj_deleteModel(m);
}

int main(int argc, char** argv) {
  const char* xml = argc > 1 ? argv[1] : "/tmp/wind_demo.xml";
  printf("Wind = 5 m/s along +x, air density 1.225, gravity off.\n");
  printf("Plate is a 0.30 x 0.30 x 0.004 m disc.\n\n");
  run(xml,  0.0,            "edge-on   (0 deg)");   // face parallel to wind -> min drag
  run(xml, 15.0*M_PI/180.0, "tilted   (15 deg)");   // tilted -> drag + lift
  run(xml, 45.0*M_PI/180.0, "tilted   (45 deg)");   // peak lift
  run(xml, 90.0*M_PI/180.0, "broadside(90 deg)");   // face into wind -> max drag
  return 0;
}
