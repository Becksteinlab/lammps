/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Gabriel Alkuino (SyracuseU)
                        Joel Clemmer (SNL)
                        Teng Zhang (SyracuseU)
------------------------------------------------------------------------- */

#include "bond_bpm_corotational.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

static constexpr double EPSILON = 1e-10;

using namespace LAMMPS_NS;
using MathConst::MY_SQRT2;
using MathConst::MY_PI;

/* ---------------------------------------------------------------------- */

static double acos_limit(double c)
{
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;
  return acos(c);
}

static double wrap_2pi(double t)
{
  if (t > MY_PI) t -= 2.0 * MY_PI;
  if (t < -MY_PI) t += 2.0 * MY_PI;
  return t;
}

/* ----------------------------------------------------------------------
   Average of two quaternions using shortest path SLERP.
------------------------------------------------------------------------- */
static void quat_average(double *qA, double *qB, double *q_out)
{
  double qB_use[4] = {qB[0], qB[1], qB[2], qB[3]};
  double dot = qA[0] * qB[0] + qA[1] * qB[1] + qA[2] * qB[2] + qA[3] * qB[3];
  if (dot < 0.0) {
    qB_use[0] = -qB[0]; qB_use[1] = -qB[1];
    qB_use[2] = -qB[2]; qB_use[3] = -qB[3];
    dot = -dot;
  }
  double s = 1.0 / sqrt(2.0 * (1.0 + dot));
  q_out[0] = (qA[0] + qB_use[0]) * s;
  q_out[1] = (qA[1] + qB_use[1]) * s;
  q_out[2] = (qA[2] + qB_use[2]) * s;
  q_out[3] = (qA[3] + qB_use[3]) * s;
  MathExtra::qnormalize(q_out);
}

/* ----------------------------------------------------------------------
   Decompose quaternion into swing and twist using Wang's formula
------------------------------------------------------------------------- */
static void decompose_swing_twist(double *quat, double *swing, double *twist)
{
  double tmp = sqrt(quat[0] * quat[0] + quat[3] * quat[3]);
  double psi;
  if (tmp < EPSILON) {
    psi = 0.0;
  } else {
    psi = 2.0 * acos_limit(quat[0] / tmp);
  }
  if (quat[3] < 0.0) psi = -psi;

  double c = cos(psi / 2.0);
  double s = sin(psi / 2.0);
  twist[0] = c; twist[1] = 0.0; twist[2] = 0.0; twist[3] = s;

  double twist_inv[4];
  MathExtra::qconjugate(twist, twist_inv);
  MathExtra::quatquat(twist_inv, quat, swing);
  MathExtra::qnormalize(swing);
}

/* ---------------------------------------------------------------------- */

BondBPMCorotational::BondBPMCorotational(LAMMPS *_lmp) :
    BondBPM(_lmp), Kr(nullptr), Ks(nullptr), Kt(nullptr), Kb(nullptr), gr(nullptr),
    gs(nullptr), gt(nullptr), gb(nullptr), Fcr(nullptr), Fcs(nullptr), Tct(nullptr),
    Tcb(nullptr)
{
  partial_flag = 1;
  smooth_flag = 1;
  normalize_flag = 0;
  writedata = 0;

  update_flag = 1;
  // History: [0]=ri_mag, [1-3]=ri_hat, [4]=gamma, [5]=theta, [6]=psi
  nhistory = 7;
  id_fix_bond_history = utils::strdup("HISTORY_BPM_COROTATIONAL");

  single_extra = 7;
  svector = new double[7];
}

/* ---------------------------------------------------------------------- */

BondBPMCorotational::~BondBPMCorotational()
{
  delete[] svector;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(Kr);
    memory->destroy(Ks);
    memory->destroy(Kt);
    memory->destroy(Kb);
    memory->destroy(Fcr);
    memory->destroy(Fcs);
    memory->destroy(Tct);
    memory->destroy(Tcb);
    memory->destroy(gr);
    memory->destroy(gs);
    memory->destroy(gt);
    memory->destroy(gb);
  }
}

/* ----------------------------------------------------------------------
  Store data for a single bond - if bond added after LAMMPS init
------------------------------------------------------------------------- */

double BondBPMCorotational::store_bond(int n, int i, int j)
{
  double delx, dely, delz, r, rinv;
  double **x = atom->x;
  tagint *tag = atom->tag;
  double **bondstore = fix_bond_history->bondstore;

  // ri points from atom with smaller tag to larger tag
  if (tag[i] < tag[j]) {
    delx = x[j][0] - x[i][0];
    dely = x[j][1] - x[i][1];
    delz = x[j][2] - x[i][2];
  } else {
    delx = x[i][0] - x[j][0];
    dely = x[i][1] - x[j][1];
    delz = x[i][2] - x[j][2];
  }

  r = sqrt(delx * delx + dely * dely + delz * delz);
  rinv = 1.0 / r;

  bondstore[n][0] = r;
  bondstore[n][1] = delx * rinv;
  bondstore[n][2] = dely * rinv;
  bondstore[n][3] = delz * rinv;
  bondstore[n][4] = 0.0;  // gamma
  bondstore[n][5] = 0.0;  // theta
  bondstore[n][6] = 0.0;  // psi

  if (i < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[i]; m++) {
      if (atom->bond_atom[i][m] == tag[j]) {
        fix_bond_history->update_atom_value(i, m, 0, r);
        fix_bond_history->update_atom_value(i, m, 1, delx * rinv);
        fix_bond_history->update_atom_value(i, m, 2, dely * rinv);
        fix_bond_history->update_atom_value(i, m, 3, delz * rinv);
        for (int a = 4; a < 7; a++)
          fix_bond_history->update_atom_value(i, m, a, 0.0);
      }
    }
  }

  if (j < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[j]; m++) {
      if (atom->bond_atom[j][m] == tag[i]) {
        fix_bond_history->update_atom_value(j, m, 0, r);
        fix_bond_history->update_atom_value(j, m, 1, delx * rinv);
        fix_bond_history->update_atom_value(j, m, 2, dely * rinv);
        fix_bond_history->update_atom_value(j, m, 3, delz * rinv);
        for (int a = 4; a < 7; a++)
          fix_bond_history->update_atom_value(j, m, a, 0.0);
      }
    }
  }

  return r;
}

/* ----------------------------------------------------------------------
  Store data for all bonds called once
------------------------------------------------------------------------- */

void BondBPMCorotational::store_data()
{
  int i, j, m, type;
  double delx, dely, delz, r, rinv;
  double **x = atom->x;
  int **bond_type = atom->bond_type;
  tagint *tag = atom->tag;

  for (i = 0; i < atom->nlocal; i++) {
    for (m = 0; m < atom->num_bond[i]; m++) {
      type = bond_type[i][m];
      if (type <= 0) continue;

      j = atom->map(atom->bond_atom[i][m]);
      if (j == -1) error->one(FLERR, "Atom missing in BPM bond");

      // ri points from smaller tag to larger tag
      if (tag[i] < tag[j]) {
        delx = x[j][0] - x[i][0];
        dely = x[j][1] - x[i][1];
        delz = x[j][2] - x[i][2];
      } else {
        delx = x[i][0] - x[j][0];
        dely = x[i][1] - x[j][1];
        delz = x[i][2] - x[j][2];
      }

      domain->minimum_image(FLERR, delx, dely, delz);
      r = sqrt(delx * delx + dely * dely + delz * delz);
      rinv = 1.0 / r;

      fix_bond_history->update_atom_value(i, m, 0, r);
      fix_bond_history->update_atom_value(i, m, 1, delx * rinv);
      fix_bond_history->update_atom_value(i, m, 2, dely * rinv);
      fix_bond_history->update_atom_value(i, m, 3, delz * rinv);
      for (int a = 4; a < 7; a++)
        fix_bond_history->update_atom_value(i, m, a, 0.0);
    }
  }

  fix_bond_history->post_neighbor();
}

double BondBPMCorotational::calculate_forces(int i1, int i2, int type,
                                                double *ri, double *rf,
                                                double *force1on2, double *torque1on2,
                                                double *torque2on1, double &ebond,
                                                double &gamma, double &theta, double &psi)
{
  double **quat = atom->quat;
  double **v = atom->v;

  double dt_inv = 1.0 / update->dt;

  double Kr_type = Kr[type];
  double Ks_type = Ks[type];

  double ri_norm = MathExtra::len3(ri);
  if (normalize_flag) {
    double ri_inv = 1.0 / ri_norm;
    Kr_type *= ri_inv;
    Ks_type *= ri_inv;
  }

  // Store prior values for damping
  double gamma_prior = gamma;
  double theta_prior = theta;
  double psi_prior = psi;

  // Current quaternions
  double q1f[4], q2f[4];
  q1f[0] = quat[i1][0]; q1f[1] = quat[i1][1]; q1f[2] = quat[i1][2]; q1f[3] = quat[i1][3];
  q2f[0] = quat[i2][0]; q2f[1] = quat[i2][1]; q2f[2] = quat[i2][2]; q2f[3] = quat[i2][3];

  // qcf: global to C frame
  double qcf[4], qcf_inv[4];
  quat_average(q1f, q2f, qcf);
  MathExtra::qconjugate(qcf, qcf_inv);

  // Radial force in C frame
  double rc[3];
  MathExtra::quatrotvec(qcf_inv, rf, rc);
  double rc_norm = MathExtra::len3(rc);
  double Fr_c[3];
  double Fr_mag = Kr_type * (rc_norm - ri_norm);

  // Damping on radial velocity
  double dv[3], dv_c[3];
  MathExtra::sub3(v[i2], v[i1], dv);
  MathExtra::quatrotvec(qcf_inv, dv, dv_c);
  double rc_hat[3];
  MathExtra::normalize3(rc, rc_hat);
  double dvdotr = MathExtra::dot3(dv_c, rc_hat);
  Fr_mag += gr[type] * dvdotr;

  MathExtra::scale3(Fr_mag, rc_hat, Fr_c);

  // Shear force in C frame
  double ri_dot_rc = MathExtra::dot3(ri, rc);
  double c = ri_dot_rc / (ri_norm * rc_norm);
  gamma = acos_limit(c);
  double dgamma = wrap_2pi(gamma - gamma_prior) * dt_inv;

  double rixrc[3], s_hat[3], t_hat[3];
  MathExtra::cross3(ri, rc, rixrc);
  double rixrc_mag = MathExtra::len3(rixrc);

  double Fst_c[3] = {0.0, 0.0, 0.0};
  double Tst_c[3] = {0.0, 0.0, 0.0};

  if (rixrc_mag > EPSILON) {
    MathExtra::normalize3(rixrc, t_hat);
    MathExtra::cross3(t_hat, rc_hat, s_hat);
    MathExtra::norm3(s_hat);

    double Fs_mag = Ks_type * rc_norm * gamma + gs[type] * rc_norm * dgamma;
    MathExtra::scale3(Fs_mag, s_hat, Fst_c);

    double Tst_mag = 0.5 * rc_norm * Fs_mag;
    MathExtra::scale3(Tst_mag, t_hat, Tst_c);
  }

  // Orientation of atom2 w.r.t. C frame
  double qu0[4];
  MathExtra::quatquat(qcf_inv, q2f, qu0);
  MathExtra::qnormalize(qu0);

  // Quaternion that rotates z to rc direction
  double qm[4];
  double temp = rc_norm + rc[2];
  if (temp < 0.0) temp = 0.0;
  qm[0] = MY_SQRT2 * 0.5 * sqrt(temp / rc_norm);

  temp = sqrt(rc[0]*rc[0] + rc[1]*rc[1]);
  if (temp > EPSILON) {
    double temp2 = rc_norm - rc[2];
    if (temp2 < 0.0) temp2 = 0.0;
    double factor = -MY_SQRT2 * 0.5 * sqrt(temp2 / rc_norm) / temp;
    qm[1] = factor * rc[1];
    qm[2] = -factor * rc[0];
  } else {
    qm[1] = 0.0;
    qm[2] = 0.0;
  }
  qm[3] = 0.0;
  if (qm[0] == 0.0 && qm[1] == 0.0 && qm[2] == 0.0) {
    qm[2] = 1.0;
  }
  MathExtra::qnormalize(qm);

  // qu = qm^-1 * qu0 * qm
  double qm_inv[4], qu[4], qtmp[4];
  MathExtra::qconjugate(qm, qm_inv);
  MathExtra::quatquat(qm_inv, qu0, qtmp);
  MathExtra::quatquat(qtmp, qm, qu);
  MathExtra::qnormalize(qu);

  // Swing-twist decomposition
  double swing2[4], twist2[4];
  decompose_swing_twist(qu, swing2, twist2);

  double qu_inv[4];
  MathExtra::qconjugate(qu, qu_inv);
  double swing1[4], twist1[4];
  decompose_swing_twist(qu_inv, swing1, twist1);

  // Combined twist = twist2 * twist1^-1
  double twist1_inv[4], twist[4];
  MathExtra::qconjugate(twist1, twist1_inv);
  MathExtra::quatquat(twist2, twist1_inv, twist);
  psi = 2.0 * acos_limit(twist[0]);

  double dpsi = wrap_2pi(psi - psi_prior) * dt_inv;

  double taxis[3];
  if (fabs(psi) > EPSILON) {
    double twist_vec[3] = {twist[1], twist[2], twist[3]};
    MathExtra::normalize3(twist_vec, taxis);
  } else {
    taxis[0] = 0.0; taxis[1] = 0.0; taxis[2] = 1.0;
  }

  // Combined swing = swing2 * swing1^-1
  double swing1_inv[4], swing[4];
  MathExtra::qconjugate(swing1, swing1_inv);
  MathExtra::quatquat(swing2, swing1_inv, swing);
  theta = 2.0 * acos_limit(swing[0]);
  double dtheta = wrap_2pi(theta - theta_prior) * dt_inv;

  double baxis[3];
  if (fabs(theta) > EPSILON) {
    double swing_vec[3] = {swing[1], swing[2], swing[3]};
    MathExtra::normalize3(swing_vec, baxis);
  } else {
    baxis[0] = 0.0; baxis[1] = 0.0; baxis[2] = 1.0;
  }

  // Torques in C' frame with damping
  double Tt_p[3], Tb_p[3];
  double Tt_mag = Kt[type] * psi + gt[type] * dpsi;
  double Tb_mag = Kb[type] * theta + gb[type] * dtheta;

  MathExtra::scale3(Tt_mag, taxis, Tt_p);
  MathExtra::scale3(Tb_mag, baxis, Tb_p);

  // Rotate torques to C frame
  double Tb_c[3], Tt_c[3];
  MathExtra::quatrotvec(qm, Tb_p, Tb_c);
  MathExtra::quatrotvec(qm, Tt_p, Tt_c);

  // Total forces and torques in C frame
  double F1_c[3], T1_c[3], T2_c[3], Ttmp[3];
  MathExtra::add3(Fr_c, Fst_c, F1_c);
  MathExtra::add3(Tb_c, Tt_c, Ttmp);
  MathExtra::add3(Tst_c, Ttmp, T1_c);
  MathExtra::sub3(Tst_c, Ttmp, T2_c);

  // Rotate to global frame
  double force2on1[3];
  MathExtra::quatrotvec(qcf, F1_c, force2on1);
  MathExtra::scale3(-1.0, force2on1, force1on2);

  MathExtra::quatrotvec(qcf, T1_c, torque2on1);
  MathExtra::quatrotvec(qcf, T2_c, torque1on2);

  // Breaking criterion (use undamped magnitudes)
  double Fs_undamped = Ks_type * rc_norm * gamma;
  double Tt_undamped = Kt[type] * fabs(psi);
  double Tb_undamped = Kb[type] * theta;

  double breaking = fabs(Fr_mag) / Fcr[type] + Fs_undamped / Fcs[type] +
                    Tb_undamped / Tcb[type] + Tt_undamped / Tct[type];
  if (breaking < 0.0) breaking = 0.0;

  // Approximate bond energy:
  ebond = 0.5 * Kr_type * (rc_norm - ri_norm) * (rc_norm - ri_norm)
	+ 0.5 * Ks_type * rc_norm * rc_norm * gamma * gamma
        + 0.5 * Kt[type] * psi * psi
        + 0.5 * Kb[type] * theta * theta;

  return breaking;
}

/* ---------------------------------------------------------------------- */

void BondBPMCorotational::compute(int eflag, int vflag)
{
  if (!fix_bond_history->stored_flag) {
    fix_bond_history->stored_flag = true;
    store_data();
  }

  if (hybrid_flag) fix_bond_history->compress_history();

  int i1, i2, itmp, n, type;
  double rf[3], ri[3];
  double ri_norm;
  double breaking, smooth;
  double force1on2[3], torque1on2[3], torque2on1[3], ebond;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  tagint *tag = atom->tag;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  double **bondstore = fix_bond_history->bondstore;

  for (n = 0; n < nbondlist; n++) {
    if (bondlist[n][2] <= 0) continue;

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    ri_norm = bondstore[n][0];

    // Order so i1 has lower tag
    if (tag[i2] < tag[i1]) {
      itmp = i1;
      i1 = i2;
      i2 = itmp;
    }

    if (ri_norm < EPSILON || std::isnan(ri_norm)) ri_norm = store_bond(n, i1, i2);

    // ri points from i1 (lower tag) to i2 (higher tag)
    ri[0] = bondstore[n][1] * ri_norm;
    ri[1] = bondstore[n][2] * ri_norm;
    ri[2] = bondstore[n][3] * ri_norm;

    // rf = x[i2] - x[i1]
    MathExtra::sub3(x[i2], x[i1], rf);

    breaking = calculate_forces(i1, i2, type, ri, rf,
                                force1on2, torque1on2, torque2on1, ebond,
                                bondstore[n][4], bondstore[n][5], bondstore[n][6]);

    if ((breaking >= 1.0) && break_flag) {
      bondlist[n][2] = 0;
      process_broken(i1, i2);
      continue;
    }

    if (smooth_flag) {
      smooth = breaking * breaking;
      smooth = 1.0 - smooth * smooth;
    } else {
      smooth = 1.0;
    }

    MathExtra::scale3(smooth, force1on2);

    if (newton_bond || i1 < nlocal) {
      f[i1][0] -= force1on2[0];
      f[i1][1] -= force1on2[1];
      f[i1][2] -= force1on2[2];

      MathExtra::scale3(smooth, torque2on1);
      torque[i1][0] += torque2on1[0];
      torque[i1][1] += torque2on1[1];
      torque[i1][2] += torque2on1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += force1on2[0];
      f[i2][1] += force1on2[1];
      f[i2][2] += force1on2[2];

      MathExtra::scale3(smooth, torque1on2);
      torque[i2][0] += torque1on2[0];
      torque[i2][1] += torque1on2[1];
      torque[i2][2] += torque1on2[2];
    }

    if (evflag)
      ev_tally_xyz(i1, i2, nlocal, newton_bond, ebond, -force1on2[0], -force1on2[1],
                   -force1on2[2], rf[0], rf[1], rf[2]);
  }

  if (hybrid_flag) fix_bond_history->uncompress_history();
}

/* ---------------------------------------------------------------------- */

void BondBPMCorotational::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(Kr, np1, "bond:Kr");
  memory->create(Ks, np1, "bond:Ks");
  memory->create(Kt, np1, "bond:Kt");
  memory->create(Kb, np1, "bond:Kb");
  memory->create(Fcr, np1, "bond:Fcr");
  memory->create(Fcs, np1, "bond:Fcs");
  memory->create(Tct, np1, "bond:Tct");
  memory->create(Tcb, np1, "bond:Tcb");
  memory->create(gr, np1, "bond:gr");
  memory->create(gs, np1, "bond:gs");
  memory->create(gt, np1, "bond:gt");
  memory->create(gb, np1, "bond:gb");

  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
   args: Kr Ks Kt Kb Fcr Fcs Tct Tcb gr gs gt gb
------------------------------------------------------------------------- */

void BondBPMCorotational::coeff(int narg, char **arg)
{
  if (narg != 13) error->all(FLERR, "Incorrect args for bond coefficients: "
                            "expected Kr Ks Kt Kb Fcr Fcs Tct Tcb gr gs gt gb");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double Kr_one = utils::numeric(FLERR, arg[1], false, lmp);
  double Ks_one = utils::numeric(FLERR, arg[2], false, lmp);
  double Kt_one = utils::numeric(FLERR, arg[3], false, lmp);
  double Kb_one = utils::numeric(FLERR, arg[4], false, lmp);
  double Fcr_one = utils::numeric(FLERR, arg[5], false, lmp);
  double Fcs_one = utils::numeric(FLERR, arg[6], false, lmp);
  double Tct_one = utils::numeric(FLERR, arg[7], false, lmp);
  double Tcb_one = utils::numeric(FLERR, arg[8], false, lmp);
  double gr_one = utils::numeric(FLERR, arg[9], false, lmp);
  double gs_one = utils::numeric(FLERR, arg[10], false, lmp);
  double gt_one = utils::numeric(FLERR, arg[11], false, lmp);
  double gb_one = utils::numeric(FLERR, arg[12], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    Kr[i] = Kr_one;
    Ks[i] = Ks_one;
    Kt[i] = Kt_one;
    Kb[i] = Kb_one;
    Fcr[i] = Fcr_one;
    Fcs[i] = Fcs_one;
    Tct[i] = Tct_one;
    Tcb[i] = Tcb_one;
    gr[i] = gr_one;
    gs[i] = gs_one;
    gt[i] = gt_one;
    gb[i] = gb_one;
    setflag[i] = 1;
    count++;

    if (Fcr[i] / Kr[i] > max_stretch) max_stretch = Fcr[i] / Kr[i];
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   check for correct settings and create fix
------------------------------------------------------------------------- */

void BondBPMCorotational::init_style()
{
  BondBPM::init_style();

  if (!atom->quat_flag || !atom->radius_flag || !atom->omega_flag)
    error->all(FLERR, "Bond bpm/corotational requires atom style bpm/sphere");
  if (comm->ghost_velocity == 0)
    error->all(FLERR, "Bond bpm/corotational requires ghost atoms store velocity");

  if (domain->dimension == 2)
    error->warning(FLERR, "Bond style bpm/corotational not intended for 2d use");
}

/* ---------------------------------------------------------------------- */

void BondBPMCorotational::settings(int narg, char **arg)
{
  BondBPM::settings(narg, arg);

  int iarg;
  for (std::size_t i = 0; i < leftover_iarg.size(); i++) {
    iarg = leftover_iarg[i];
    if (strcmp(arg[iarg], "smooth") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal bond bpm command, missing option for smooth");
      smooth_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      i += 1;
    } else if (strcmp(arg[iarg], "normalize") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal bond bpm command, missing option for normalize");
      normalize_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      i += 1;
    } else {
      error->all(FLERR, "Illegal bond bpm command, invalid argument {}", arg[iarg]);
    }
  }

  if (smooth_flag && !break_flag)
    error->all(FLERR, "Illegal bond bpm command, must turn off smoothing with break no option");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondBPMCorotational::write_restart(FILE *fp)
{
  BondBPM::write_restart(fp);
  write_restart_settings(fp);

  fwrite(&Kr[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&Ks[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&Kt[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&Kb[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&Fcr[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&Fcs[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&Tct[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&Tcb[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gr[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gs[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gt[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gb[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondBPMCorotational::read_restart(FILE *fp)
{
  BondBPM::read_restart(fp);
  read_restart_settings(fp);
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &Kr[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &Ks[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &Kt[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &Kb[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &Fcr[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &Fcs[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &Tct[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &Tcb[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gr[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gs[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gt[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gb[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&Kr[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Ks[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Kt[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Kb[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Fcr[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Fcs[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Tct[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Tcb[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gr[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gs[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gt[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gb[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void BondBPMCorotational::write_restart_settings(FILE *fp)
{
  fwrite(&smooth_flag, sizeof(int), 1, fp);
  fwrite(&normalize_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void BondBPMCorotational::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &smooth_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &normalize_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&smooth_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&normalize_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double BondBPMCorotational::single(int type, double rsq, int i, int j, double &fforce)
{
  if (type <= 0) return 0.0;

  int flipped = 0;
  if (atom->tag[j] < atom->tag[i]) {
    int itmp = i;
    i = j;
    j = itmp;
    flipped = 1;
  }

  double ri_norm, gamma, theta, psi;
  double ri[3], rf[3];
  for (int n = 0; n < atom->num_bond[i]; n++) {
    if (atom->bond_atom[i][n] == atom->tag[j]) {
      ri_norm = fix_bond_history->get_atom_value(i, n, 0);
      ri[0] = fix_bond_history->get_atom_value(i, n, 1) * ri_norm;
      ri[1] = fix_bond_history->get_atom_value(i, n, 2) * ri_norm;
      ri[2] = fix_bond_history->get_atom_value(i, n, 3) * ri_norm;
      gamma = fix_bond_history->get_atom_value(i, n, 4);
      theta = fix_bond_history->get_atom_value(i, n, 5);
      psi = fix_bond_history->get_atom_value(i, n, 6);
    }
  }

  double **x = atom->x;
  MathExtra::sub3(x[j], x[i], rf);

  double force1on2[3], torque1on2[3], torque2on1[3], ebond;
  double breaking = calculate_forces(i, j, type, ri, rf,
                                     force1on2, torque1on2, torque2on1, ebond,
                                     gamma, theta, psi);

  double rf_hat[3];
  double rf_mag = sqrt(rsq);
  MathExtra::scale3(1.0/rf_mag, rf, rf_hat);
  fforce = MathExtra::dot3(force1on2, rf_hat);

  double smooth = 1.0;
  if (smooth_flag) {
    smooth = breaking * breaking;
    smooth = 1.0 - smooth * smooth;
    fforce *= smooth;
  }

  MathExtra::scale3(smooth, force1on2);
  svector[0] = ri_norm;
  if (flipped) {
    svector[1] = -ri[0];
    svector[2] = -ri[1];
    svector[3] = -ri[2];
    svector[4] = -force1on2[0];
    svector[5] = -force1on2[1];
    svector[6] = -force1on2[2];
  } else {
    svector[1] = ri[0];
    svector[2] = ri[1];
    svector[3] = ri[2];
    svector[4] = force1on2[0];
    svector[5] = force1on2[1];
    svector[6] = force1on2[2];
  }

  return ebond;
}
