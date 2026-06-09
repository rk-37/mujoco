"""Interceptor simulation harness (Python).

Drives the scenario and feeds the drone's onboard camera (720p) + IMU/altitude
sensors to controller.control(), which returns the four motor thrusts.

Pick the control strategy with --controller NAME (a module in controllers/).

Usage:
    .venv/bin/python   run_intercept.py [--controller hover] [--record out.mp4]
    .venv/bin/mjpython run_intercept.py --live [--controller hover]
"""

import sys
import importlib
import numpy as np
import mujoco

XML = "/Users/rihards/src/mujoco/experiments/interceptor/intercept.xml"


def load_controller(name):
    """Import controllers/<name>.py and return its control() function."""
    mod = importlib.import_module(f"controllers.{name}")
    print(f"[controller] using '{name}'")
    return mod.control

# scenario constants (mirror the C++ harness)
V_TGT, ALT, X0, XEND, YOFF = 15.0, 30.0, -150.0, 150.0, 30.0
OBS = np.array([0.0, 0.0, 1.7])          # ground observer eye
LAUNCH_FROM = np.array([2.0, 0.0, 1.2])  # throw origin (beside observer)
THROW_SPEED = 13.41                      # 30 mph
LAUNCH_RANGE = 80.0                      # throw when target within this of the observer
AIM_NOISE = np.radians(4.0)             # directional noise on the throw
CONTROL_HZ = 30
FPV_W, FPV_H = 1280, 720                  # 720p feed to the controller
MAIN_FOVY = 18.0                          # binocular-ish zoom for the main view


class Sim:
    def __init__(self, control_fn):
        self.control_fn = control_fn
        self.m = mujoco.MjModel.from_xml_path(XML)
        self.d = mujoco.MjData(self.m)
        nid = lambda t, n: mujoco.mj_name2id(self.m, t, n)
        OBJ = mujoco.mjtObj
        self.t_bid = nid(OBJ.mjOBJ_BODY, "target")
        self.i_bid = nid(OBJ.mjOBJ_BODY, "interceptor")
        self.t_q = self.m.jnt_qposadr[self.m.body_jntadr[self.t_bid]]
        self.t_v = self.m.body_dofadr[self.t_bid]
        self.i_q = self.m.jnt_qposadr[self.m.body_jntadr[self.i_bid]]
        self.i_v = self.m.body_dofadr[self.i_bid]
        self.t_act = [nid(OBJ.mjOBJ_ACTUATOR, n) for n in ("m_px", "m_nx", "m_py", "m_ny")]
        self.i_act = [nid(OBJ.mjOBJ_ACTUATOR, n) for n in ("i_m_px", "i_m_nx", "i_m_py", "i_m_ny")]
        self.fpv = nid(OBJ.mjOBJ_CAMERA, "fpv")
        self.tmass = self.m.body_mass[self.t_bid]
        self.launched = False
        self.theta = 0.0
        self.Tper = 0.0
        self.rng = np.random.default_rng()
        self._trim()
        self.reset()

    # ---- target trim (bisection on pitch, Newton on thrust) ----
    def _eval(self, theta, mg):
        T = mg / np.cos(theta)
        self.d.qpos[self.t_q:self.t_q + 3] = (0, 0, ALT)
        self.d.qpos[self.t_q + 3:self.t_q + 7] = (np.cos(theta / 2), 0, np.sin(theta / 2), 0)
        self.d.qvel[self.t_v:self.t_v + 6] = (V_TGT, 0, 0, 0, 0, 0)
        for _ in range(3):
            for a in self.t_act:
                self.d.ctrl[a] = T / 4
            mujoco.mj_forward(self.m, self.d)
            az = self.d.qacc[self.t_v + 2]
            if abs(az) < 1e-12:
                break
            T -= az * self.tmass / np.cos(theta)
        return self.d.qacc[self.t_v + 0], T / 4

    def _trim(self):
        mg = self.tmass * (-self.m.opt.gravity[2])
        lo, hi = np.radians(-25), np.radians(25)
        axlo, _ = self._eval(lo, mg)
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            ax, Tper = self._eval(mid, mg)
            if (ax > 0) == (axlo > 0):
                lo, axlo = mid, ax
            else:
                hi = mid
        self.theta, self.Tper = 0.5 * (lo + hi), self._eval(0.5 * (lo + hi), mg)[1]
        print(f"[trim] target theta={np.degrees(self.theta):.2f} deg  T/rotor={self.Tper:.3f} N")

    # ---- state setters ----
    def _set_target_start(self):
        self.d.qpos[self.t_q:self.t_q + 3] = (X0, YOFF, ALT)
        self.d.qpos[self.t_q + 3:self.t_q + 7] = (np.cos(self.theta / 2), 0, np.sin(self.theta / 2), 0)
        self.d.qvel[self.t_v:self.t_v + 6] = (V_TGT, 0, 0, 0, 0, 0)

    def _stow(self):
        self.launched = False
        self.d.qpos[self.i_q:self.i_q + 3] = LAUNCH_FROM
        self.d.qpos[self.i_q + 3:self.i_q + 7] = (1, 0, 0, 0)
        self.d.qvel[self.i_v:self.i_v + 6] = 0
        for a in self.i_act:
            self.d.ctrl[a] = 0

    def reset(self):
        mujoco.mj_resetData(self.m, self.d)
        self._set_target_start()
        self._stow()
        mujoco.mj_forward(self.m, self.d)

    def _launch(self):
        tp = self.d.xpos[self.t_bid].copy()
        dir = tp - LAUNCH_FROM
        dir /= np.linalg.norm(dir)
        dir = dir + AIM_NOISE * self.rng.standard_normal(3)
        dir /= np.linalg.norm(dir)
        self.aim = dir
        # level body: motors + up-camera point straight up, hover thrust cancels gravity
        self.d.qpos[self.i_q + 3:self.i_q + 7] = (1, 0, 0, 0)
        self.d.qvel[self.i_v:self.i_v + 3] = THROW_SPEED * dir
        self.d.qvel[self.i_v + 3:self.i_v + 6] = 0
        self.launched = True
        print(f"[launch] dir={np.round(dir,2)} toward target {np.round(tp,0)}")

    def sensors(self):
        def s(name):
            i = mujoco.mj_name2id(self.m, mujoco.mjtObj.mjOBJ_SENSOR, name)
            a, n = self.m.sensor_adr[i], self.m.sensor_dim[i]
            return np.array(self.d.sensordata[a:a + n])
        return {"accel": s("i_acc"), "gyro": s("i_gyro"), "alt": float(s("i_pos")[2])}

    def step_control(self, frame):
        """one control tick: target trim + interceptor control() + physics steps."""
        for a in self.t_act:
            self.d.ctrl[a] = self.Tper
        if not self.launched:
            self._stow()
            if np.linalg.norm(self.d.xpos[self.t_bid] - LAUNCH_FROM) < LAUNCH_RANGE:
                self._launch()
        else:
            cmd = self.control_fn(frame, self.sensors())
            for k, a in enumerate(self.i_act):
                self.d.ctrl[a] = float(cmd[k])
        nstep = max(1, round((1.0 / CONTROL_HZ) / self.m.opt.timestep))
        for _ in range(nstep):
            mujoco.mj_step(self.m, self.d)
        if self.d.qpos[self.t_q] >= XEND:
            self.reset()
            return True   # completed a pass
        return False


# ----- camera-facing dotted square markers (for the main view) -----
def add_marker(scene, c, half, rgb):
    f = scene.camera[0].forward
    u = scene.camera[0].up
    R = np.cross(f, u); R /= (np.linalg.norm(R) + 1e-9)
    U = np.array(u) / (np.linalg.norm(u) + 1e-9)
    cor = [c - half*R + half*U, c + half*R + half*U, c + half*R - half*U, c - half*R - half*U]
    N = 9
    rad = max(0.05, half * 0.05)
    mat = np.eye(3).flatten()
    rgba = np.array([rgb[0], rgb[1], rgb[2], 0.55], dtype=np.float32)
    for e in range(4):
        A, B = cor[e], cor[(e + 1) % 4]
        for i in range(N):
            if scene.ngeom >= scene.maxgeom:
                return
            p = A + (B - A) * (i / N)
            g = scene.geoms[scene.ngeom]
            mujoco.mjv_initGeom(g, mujoco.mjtGeom.mjGEOM_SPHERE,
                                np.array([rad, rad, rad]), p.astype(float), mat, rgba)
            g.category = mujoco.mjtCatBit.mjCAT_DECOR
            scene.ngeom += 1


def resize_nn(img, w, h):
    ih, iw = img.shape[:2]
    ys = (np.arange(h) * ih // h)
    xs = (np.arange(w) * iw // w)
    return img[ys][:, xs]


def main_camera(sim):
    cam = mujoco.MjvCamera()
    cam.type = mujoco.mjtCamera.mjCAMERA_FREE
    return cam


def aim_main_camera(sim, cam):
    tp = sim.d.xpos[sim.t_bid]
    vv = tp - OBS
    rng = np.linalg.norm(vv)
    cam.lookat[:] = tp
    cam.distance = rng
    cam.azimuth = np.degrees(np.arctan2(vv[1], vv[0]))
    cam.elevation = np.degrees(np.arcsin(vv[2] / rng))
    return rng


def record(path, control_fn, seconds=22.0):
    import imageio
    sim = Sim(control_fn)
    sim.m.vis.global_.fovy = MAIN_FOVY
    rend = mujoco.Renderer(sim.m, FPV_H, FPV_W)   # one renderer, reused for both views
    cam = main_camera(sim)
    writer = imageio.get_writer(path, fps=CONTROL_HZ, macro_block_size=None)
    nframes = int(seconds * CONTROL_HZ)
    pw, ph, mgn = 320, 180, 12
    for _ in range(nframes):
        # FPV (720p) -> controller input
        rend.update_scene(sim.d, camera=sim.fpv)
        fpv = rend.render().copy()
        # main view (observer, binocular zoom, dotted markers)
        aim_main_camera(sim, cam)
        rend.update_scene(sim.d, camera=cam)
        add_marker(rend.scene, sim.d.xpos[sim.t_bid].copy(), 1.6, (1, 0, 0))
        add_marker(rend.scene, sim.d.xpos[sim.i_bid].copy(), 0.7, (0, 1, 0))
        main = rend.render().copy()
        # picture-in-picture of the FPV feed, bottom-right with white border
        small = resize_nn(fpv, pw, ph)
        r0, c0 = FPV_H - ph - mgn, FPV_W - pw - mgn
        main[r0 - 3:r0 + ph + 3, c0 - 3:c0 + pw + 3] = 255
        main[r0:r0 + ph, c0:c0 + pw] = small
        writer.append_data(main)
        sim.step_control(fpv)
    writer.close()
    print(f"[record] wrote {path}  ({nframes} frames @ {CONTROL_HZ} fps)")


def live(control_fn):
    import mujoco.viewer
    sim = Sim(control_fn)
    sim.m.vis.global_.fovy = MAIN_FOVY
    rend = mujoco.Renderer(sim.m, FPV_H, FPV_W)
    with mujoco.viewer.launch_passive(sim.m, sim.d) as v:
        v.cam.type = mujoco.mjtCamera.mjCAMERA_FREE
        while v.is_running():
            rend.update_scene(sim.d, camera=sim.fpv)
            fpv = rend.render().copy()
            sim.step_control(fpv)
            v.sync()


def _arg(flag, default):
    if flag in sys.argv:
        i = sys.argv.index(flag)
        if i + 1 < len(sys.argv):
            return sys.argv[i + 1]
    return default


if __name__ == "__main__":
    name = _arg("--controller", "hover")
    ctrl = load_controller(name)
    if "--live" in sys.argv:
        live(ctrl)
    else:
        record(_arg("--record", "/tmp/intercept.mp4"), ctrl)
