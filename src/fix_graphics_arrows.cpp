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

#include "fix_graphics_arrows.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "dump_image.h"
#include "error.h"
#include "group.h"
#include "input.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { NONE, DIPOLE, FORCE, VELOCITY, VARIABLE, CHUNK };

/* ---------------------------------------------------------------------- */

FixGraphicsArrows::FixGraphicsArrows(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), val{0.0, 0.0, 0.0}, xstr(nullptr), ystr(nullptr), zstr(nullptr),
    id_chunk(nullptr), id_pos(nullptr), id_vec(nullptr), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix graphics/arrows", error);

  // parse mandatory arg

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, "Illegal fix graphics/arrows nevery value");
  global_freq = nevery;
  dynamic_group_allow = 1;

  mode = NONE;
  varflag = 0;

  if (strcmp(arg[4], "dipole") == 0) {
    mode = DIPOLE;
  } else if (strcmp(arg[4], "force") == 0) {
    mode = FORCE;
  } else if (strcmp(arg[4], "velocity") == 0) {
    mode = VELOCITY;
  } else if (strcmp(arg[4], "variable") == 0) {
    mode = VARIABLE;
    if (narg < 10) utils::missing_cmd_args(FLERR, "fix graphics/arrows variable", error);
    if (strstr(arg[5], "v_") == arg[5]) {
      varflag = 1;
      xstr = utils::strdup(arg[5] + 2);
    } else {
      val[0] = utils::numeric(FLERR, arg[5], false, lmp);
    }
    if (strstr(arg[6], "v_") == arg[6]) {
      varflag = 1;
      ystr = utils::strdup(arg[6] + 2);
    } else {
      val[1] = utils::numeric(FLERR, arg[6], false, lmp);
    }
    if (strstr(arg[7], "v_") == arg[7]) {
      varflag = 1;
      ystr = utils::strdup(arg[7] + 2);
    } else {
      val[2] = utils::numeric(FLERR, arg[7], false, lmp);
    }
    radius = utils::numeric(FLERR, arg[8], false, lmp);
    if (radius <= 0.0) error->all(FLERR, 6, "Arrow radius must be > 0");
  } else if (strcmp(arg[4], "chunk") == 0) {
    mode = CHUNK;
    if (narg < 11) utils::missing_cmd_args(FLERR, "fix graphics/arrows chunk", error);

    scale = utils::numeric(FLERR, arg[8], false, lmp);
    radius = utils::numeric(FLERR, arg[9], false, lmp);
    if (radius <= 0.0) error->all(FLERR, 9, "Arrow radius must be > 0");
  } else {
    error->all(FLERR, 4, "Unknown fix graphics/arrows keyword: {}", arg[4]);
  }

  // we have the same arguments for these modes
  if ((mode == DIPOLE) || (mode == FORCE) || (mode == VELOCITY)) {
    scale = utils::numeric(FLERR, arg[5], false, lmp);
    radius = utils::numeric(FLERR, arg[6], false, lmp);
    if (radius <= 0.0) error->all(FLERR, 6, "Arrow radius must be > 0");
  }

  if ((mode == DIPOLE) && !atom->mu_flag)
    error->all(FLERR, 4, "Fix graphics/arrows dipole mode requires atom attribute mu");

  // checks
  numobjs = 0;
}

/* ---------------------------------------------------------------------- */

FixGraphicsArrows::~FixGraphicsArrows()
{
  memory->destroy(imgobjs);
  memory->destroy(imgparms);

  delete[] xstr;
  delete[] ystr;
  delete[] zstr;
  delete[] id_chunk;
  delete[] id_pos;
  delete[] id_vec;
}

/* ---------------------------------------------------------------------- */

int FixGraphicsArrows::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsArrows::init()
{
  if (xstr) {
    int ivar = input->variable->find(xstr);
    if (ivar < 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "Variable name {} for fix graphics/arrows x valaue does not exist", xstr);
    if (input->variable->atomstyle(ivar) == 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix graphics/arrows variable {} is not atom-style variable", xstr);
    xvar = ivar;
  }
  if (ystr) {
    int ivar = input->variable->find(ystr);
    if (ivar < 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "Variable name {} for fix graphics/arrows y valaue does not exist", ystr);
    if (input->variable->atomstyle(ivar) == 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix graphics/arrows variable {} is not atom-style variable", ystr);
    yvar = ivar;
  }
  if (zstr) {
    int ivar = input->variable->find(zstr);
    if (ivar < 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "Variable name {} for fix graphics/arrows z valaue does not exist", zstr);
    if (input->variable->atomstyle(ivar) == 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix graphics/arrows variable {} is not atom-style variable", zstr);
    zvar = ivar;
  }
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGraphicsArrows::end_of_step()
{
  memory->destroy(imgobjs);
  memory->destroy(imgparms);

  // evaluate variable if necessary, wrap with clear/add

  if (varflag) modify->clearstep_compute();

  const auto *const *const x = atom->x;
  const auto *const *const f = atom->f;
  const auto *const *const v = atom->v;
  const auto *const *const mu = atom->mu;
  const auto *const mask = atom->mask;
  const auto *const type = atom->type;
  const auto nlocal = atom->nlocal;

  // count (local) number of arrows

  if (mode != CHUNK) {
    int n = 0;
    for (int i = 0; i < nlocal; ++i)
      if (mask[i] & groupbit) ++n;

    numobjs = n;
    memory->create(imgobjs, numobjs, "fix_graphics_arrows:imgobjs");
    memory->create(imgparms, numobjs, 9, "fix_graphics_arrows:imgparms");

    n = 0;
    for (int i = 0; i < nlocal; ++i)
      if (mask[i] & groupbit) {
        imgobjs[n] = DumpImage::ARROW;
        imgparms[n][0] = type[i];
        imgparms[n][1] = x[i][0];
        imgparms[n][2] = x[i][1];
        imgparms[n][3] = x[i][2];
        if (mode == DIPOLE) {
          imgparms[n][4] = mu[i][0];
          imgparms[n][5] = mu[i][1];
          imgparms[n][6] = mu[i][2];
        } else if (mode == FORCE) {
          imgparms[n][4] = f[i][0];
          imgparms[n][5] = f[i][1];
          imgparms[n][6] = f[i][2];
        } else if (mode == VELOCITY) {
          imgparms[n][4] = v[i][0];
          imgparms[n][5] = v[i][1];
          imgparms[n][6] = v[i][2];
        }
        imgparms[n][7] = scale;
        imgparms[n][8] = 2.0 * radius;
        ++n;
      }
  }

  if (varflag) modify->addstep_compute((update->ntimestep / nevery) * nevery + nevery);
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image on universe root only
------------------------------------------------------------------------- */

int FixGraphicsArrows::image(int *&objs, double **&parms)
{
  objs = imgobjs;
  parms = imgparms;
  return numobjs;
}
