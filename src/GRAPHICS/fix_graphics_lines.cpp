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

#include "fix_graphics_lines.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "graphics.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGraphicsLines::FixGraphicsLines(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), id_cprop(nullptr), cprop(nullptr), id_fave(nullptr), fave(nullptr),
    id_fstore(nullptr), fstore(nullptr), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix graphics/lines", error);

  // parse mandatory args

  nave = utils::inumeric(FLERR, arg[3], false, lmp);
  nrepeat = utils::inumeric(FLERR, arg[4], false, lmp);
  nevery = utils::inumeric(FLERR, arg[5], false, lmp);
  nmax = utils::inumeric(FLERR, arg[6], false, lmp);

  // fix settings
  scalar_flag = 1;
  extscalar = 0;
  time_depend = 1;
  dynamic_group_allow = 0;    // there is no clean way to use dynamic groups for this

  nvalues = ivalue = numobjs = 0;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsLines::post_constructor()
{
  const auto *gname = group->names[igroup];
  id_cprop = utils::strdup(id + std::string("_COMPUTE_PROPERTY_ATOM"));
  cprop = modify->add_compute(fmt::format("{0} {1} property/atom xu yu zu", id_cprop, gname));
  id_fave = utils::strdup(id + std::string("_FIX_AVE_ATOM"));
  fave = modify->add_fix(fmt::format("{0} {1} ave/atom {2} {3} {4} c_{5}[1] c_{5}[2] c_{5}[3]",
                                     id_fave, gname, nave, nrepeat, nevery, id_cprop));
  id_fstore = utils::strdup(id + std::string("_FIX_STORE_ATOM"));
  fstore = modify->add_fix(fmt::format("{0} {1} STORE/ATOM {2} 0 0 1", id_fstore, gname, nmax));

  // is there a better way to do the following instead of this ugly hack?
  // swap our position with fix ave/atom, so it is executed first and we can access its data
  // we don't need to swap fmask since both fixes are END_OF_STEP only ...
  int ifave = modify->find_fix(id_fave);
  int ilines = modify->find_fix(id);
  modify->fix[ifave] = this;
  modify->fix[ilines] = fave;

  // ... but we also need to change the fix index for fix ave/atom in the Atom::grow() callback list
  for (int i = 0; i < atom->nextra_grow; ++i)
    if (atom->extra_grow[i] == ifave) atom->extra_grow[i] = ilines;
}

/* ---------------------------------------------------------------------- */

FixGraphicsLines::~FixGraphicsLines()
{
  if (modify->ncompute && modify->get_compute_by_id(id_cprop)) modify->delete_compute(id_cprop);
  if (modify->nfix && modify->get_fix_by_id(id_fave)) modify->delete_fix(id_fave);
  if (modify->nfix && modify->get_fix_by_id(id_fstore)) modify->delete_fix(id_fstore);
  delete[] id_cprop;
  delete[] id_fave;
  delete[] id_fstore;
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ---------------------------------------------------------------------- */

int FixGraphicsLines::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsLines::end_of_step()
{
  memory->destroy(imgobjs);
  memory->destroy(imgparms);

  auto num = group->count(igroup);
  numobjs = num * nmax;

#if 0
  for (int n = 0; n < numobjs; ++n) {
    imgobjs[n] = Graphics::ARROW;
    imgparms[n][0] = 1;
    imgparms[n][1] = pdata[n][0];
    imgparms[n][2] = pdata[n][1];
    imgparms[n][3] = pdata[n][2];
    imgparms[n][4] = vdata[n][0];
    imgparms[n][5] = vdata[n][1];
    imgparms[n][6] = vdata[n][2];
    imgparms[n][7] = scale;
    imgparms[n][8] = 2.0 * radius;
    imgparms[n][9] = 0.2;
  }

  if (varflag || (mode == CHUNK))
    modify->addstep_compute((update->ntimestep / nevery) * nevery + nevery);
#endif
}

/* ----------------------------------------------------------------------
   current length of lines
------------------------------------------------------------------------- */

double FixGraphicsLines::compute_scalar()
{
  return nvalues;
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image
------------------------------------------------------------------------- */

int FixGraphicsLines::image(int *&objs, double **&parms)
{
  objs = imgobjs;
  parms = imgparms;
  return numobjs;
}
