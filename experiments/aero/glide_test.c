#include <stdio.h>
#include <mujoco/mujoco.h>

int main(void) {
  char err[1000] = "";
  mjModel* m = mj_loadXML("/tmp/glider.xml", NULL, err, sizeof(err));
  if (!m) { printf("load error: %s\n", err); return 1; }
  mjData* d = mj_makeData(m);

  // load the launch keyframe (initial pose + 14 m/s forward velocity)
  mj_resetDataKeyframe(m, d, 0);

  printf(" time     x(m)    z(m)    speed   pitch(deg)\n");
  double next = 0;
  while (d->time < 6.0) {
    if (d->time >= next) {
      double vx = d->qvel[0], vz = d->qvel[2];
      double speed = mju_sqrt(vx*vx + vz*vz);
      // pitch from body quaternion (rotation about y)
      double pitch = 2.0 * mju_atan2(d->qpos[5], d->qpos[3]) * 180.0 / 3.14159265;
      printf(" %4.2f  %7.2f %7.2f  %6.2f   %+7.1f\n",
             d->time, d->qpos[0], d->qpos[2], speed, pitch);
      next += 0.5;
    }
    mj_step(m, d);
    if (d->qpos[2] < 0.05) { printf(" --- landed ---\n"); break; }
  }
  printf("\nGlide distance: %.1f m from %.1f m drop.\n", d->qpos[0], 3.0 - d->qpos[2]);
  mj_deleteData(d);
  mj_deleteModel(m);
  return 0;
}
