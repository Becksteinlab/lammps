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

#include "fix_graphics_replica.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "dump_image.h"
#include "error.h"
#include "input.h"
#include "lattice.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>
#include <utility>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGraphicsReplica::FixGraphicsReplica(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix graphics/replica", error);

  // parse mandatory arg

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix graphics/replica nevery value {}", nevery);
  global_freq = nevery;
  dynamic_group_allow = 1;

  numobjs = 0;

  int iarg = 4;
  while (iarg < narg) {
    ++iarg;
  }
  
  memory->create(imgobjs, numobjs, "fix_graphics:imgobjs");
  memory->create(imgparms, numobjs, 8, "fix_graphics:imgparms");
}

/* ---------------------------------------------------------------------- */

FixGraphicsReplica::~FixGraphicsReplica()
{
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ---------------------------------------------------------------------- */

int FixGraphicsReplica::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsReplica::init()
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGraphicsReplica::end_of_step()
{
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image
------------------------------------------------------------------------- */

int FixGraphicsReplica::image(int *&objs, double **&parms)
{
  objs = imgobjs;
  parms = imgparms;
  return numobjs;
}
