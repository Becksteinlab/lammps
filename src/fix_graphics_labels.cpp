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

#include "fix_graphics_labels.h"

#include "atom.h"
#include "comm.h"
#include "dump_image.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
#include "text_file_reader.h"
#include "tokenizer.h"
#include "update.h"
#include "variable.h"

#include <cstring>

#ifdef LAMMPS_JPEG
#include <jpeglib.h>
#endif

#ifdef LAMMPS_PNG
#include <csetjmp>
#include <png.h>
#include <zlib.h>
#endif

using namespace LAMMPS_NS;
using namespace FixConst;

namespace {

// read image into buffer that is locally allocated with new
// return null pointer if incompatible format or not supported

unsigned char *read_image(FILE *fp, int &width, int &height)
{
  if (!fp) return nullptr;
  unsigned char *pixmap = nullptr;

#if defined(LAMMPS_JPEG)
#else
  return nullptr;
#endif
#if defined(LAMMPS_PNG)
#else
  return nullptr;
#endif

  // read file in NetPBM binary format
  char buffer[128];
  char *ptr = fgets(buffer, 128, fp);
  if (!ptr || (strlen(buffer) < 3) || (buffer[0] != 'P')) return nullptr;

  // detect binary versus ASCII variant
  bool binary = true;
  if (buffer[1] == '3')
    binary = false;
  else if (buffer[1] != '6')
    return nullptr;

  // skip over optional comments
  do {
    char *ptr = fgets(buffer, 128, fp);
  } while (buffer[0] == '#');

  int rv = sscanf(buffer, "%d%d", &width, &height);
  if (rv != 2) return nullptr;

  int tmp = 0;
  ptr = fgets(buffer, 128, fp);
  rv = sscanf(buffer, "%d", &tmp);
  if ((rv != 1) || (tmp != 255)) return nullptr;

  pixmap = new unsigned char[3 * width * height];
  if (binary) {
    // read raw data directly into buffer in the expected order of lines
    // this is the inverse of what Image::write_PPM() does
    for (int y = height - 1; y >= 0; --y) {
      rv = fread(&pixmap[y * width * 3], 3, width, fp);
      if (rv != width) {
        delete[] pixmap;
        return nullptr;
      }
    }
  } else {
    // read file line-by-line and store three RGB values at a time
    auto reader = TextFileReader(fp, "NetPBM ASCII pixmap");
    try {
      int y = height -1;
      int x = 0;
      auto *line = reader.next_line();
      while (line) {
        auto values = ValueTokenizer(line);

        while (values.has_next()) {
          pixmap[y * width * 3 + 3 * x] = values.next_int();
          pixmap[y * width * 3 + 3 * x + 1] = values.next_int();
          pixmap[y * width * 3 + 3 * x + 2] = values.next_int();
          ++x;
          // next line of pixels
          if (x >= width) {
            --y;
            x = 0;
          }
        }
        line = reader.next_line();
      }
    } catch (std::exception &e) {
      delete[] pixmap;
      return nullptr;
    }
  }

  return pixmap;
}
}    // namespace

/* ---------------------------------------------------------------------- */

FixGraphicsLabels::FixGraphicsLabels(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix graphics/labels", error);

  // parse mandatory arg

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix graphics/labels nevery value");
  global_freq = nevery;
  dynamic_group_allow = 1;

  // set defaults
  numobjs = 0;
  varflag = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "image") == 0) {
      if (iarg + 5 > narg) utils::missing_cmd_args(FLERR, "fix graphics/labels", error);

      // clang-format off
      PixmapInfo pix{{0.0, 0.0, 0.0}, 0, 0, nullptr, {-1.0, -1.0, -1.0}, 1.0,
        -1, -1, -1, -1, nullptr, nullptr, nullptr, nullptr};
      // clang-format on

      // read and store image file with pixmap only on MPI rank 0
      if (comm->me == 0) {
        FILE *fp = fopen(arg[iarg + 1], "rb");
        if (!fp)
          error->one(FLERR, iarg + 1, "Cannot open fix graphics/labels image file {}: {}",
                     arg[iarg + 1], utils::getsyserror());
        pix.pixmap = read_image(fp, pix.width, pix.height);
        if (!pix.pixmap)
          error->one(FLERR, iarg + 1,
                     "Reading open fix graphics/labels image file {} failed.\n"
                     "                Unsupported file format or broken file",
                     arg[iarg + 1]);
      }

      if (strstr(arg[iarg + 2], "v_") == arg[iarg + 2]) {
        varflag = 1;
        pix.xstr = utils::strdup(arg[iarg + 2] + 2);
      } else
        pix.pos[0] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      if (strstr(arg[iarg + 3], "v_") == arg[iarg + 3]) {
        varflag = 1;
        pix.ystr = utils::strdup(arg[iarg + 3] + 2);
      } else
        pix.pos[1] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      if (strstr(arg[iarg + 4], "v_") == arg[iarg + 4]) {
        varflag = 1;
        pix.zstr = utils::strdup(arg[iarg + 4] + 2);
      } else
        pix.pos[2] = utils::numeric(FLERR, arg[iarg + 4], false, lmp);

      iarg += 5;

      // check remaining arguments for optional image arguments
      while (iarg < narg) {
        // next argument is next keyword; exit loop
        if ((strcmp(arg[iarg], "image") == 0) || (strcmp(arg[iarg], "text") == 0)) break;

        if (strcmp(arg[iarg], "scale") == 0) {
          if (iarg + 2 > narg)
            utils::missing_cmd_args(FLERR, "fix graphics/labels image scale", error);
          if (strstr(arg[iarg + 1], "v_") == arg[iarg + 1]) {
            varflag = 1;
            pix.sstr = utils::strdup(arg[iarg + 1] + 2);
          } else
            pix.scale = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
          if (pix.scale <= 0.0)
            error->all(FLERR, iarg + 1, "Invalid fix graphics/labels image scale value: {}",
                       pix.scale);
          iarg += 2;
        } else if (strcmp(arg[iarg], "transcolor") == 0) {
          if (iarg + 2 > narg)
            utils::missing_cmd_args(FLERR, "fix graphics/labels image transcolor", error);
          if (strcmp(arg[iarg + 1], "auto") == 0) {
            pix.transcolor[0] = pix.pixmap[0] / 255.0;
            pix.transcolor[1] = pix.pixmap[1] / 255.0;
            pix.transcolor[2] = pix.pixmap[2] / 255.0;
          } else if (strcmp(arg[iarg + 1], "none") == 0) {
            pix.transcolor[0] = -1.0;
            pix.transcolor[1] = -1.0;
            pix.transcolor[2] = -1.0;
          } else {
            auto rgb = ValueTokenizer(arg[iarg + 1], "/");
            try {
              pix.transcolor[0] = rgb.next_int() / 255.0;
              pix.transcolor[1] = rgb.next_int() / 255.0;
              pix.transcolor[2] = rgb.next_int() / 255.0;
              if (rgb.has_next()) throw TokenizerException("Extra token", rgb.next_string());
            } catch (TokenizerException &e) {
              error->all(FLERR, iarg + 1, "Error parsing RGB color value {}: {}", arg[iarg + 1],
                         e.what());
            }
          }
          iarg += 2;
        } else {
          error->all(FLERR, iarg, "Unknown fix graphics/labels image keyword: {}", arg[iarg]);
        }
      }
      pixmaps.emplace_back(pix);
    } else if (strcmp(arg[iarg], "text") == 0) {
      iarg += 1;
    } else {
      error->all(FLERR, iarg, "Unknown fix graphics/labels keyword: {}", arg[iarg]);
    }
  }
}

/* ---------------------------------------------------------------------- */

FixGraphicsLabels::~FixGraphicsLabels()
{
  for (auto &pix : pixmaps) {
    delete[] pix.pixmap;
    delete[] pix.xstr;
    delete[] pix.ystr;
    delete[] pix.zstr;
  }

  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ---------------------------------------------------------------------- */

int FixGraphicsLabels::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsLabels::init()
{
  for (auto &pix : pixmaps) {
    if (pix.xstr) {
      int ivar = input->variable->find(pix.xstr);
      if (ivar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix graphics/labels does not exist", pix.xstr);
      if (input->variable->equalstyle(ivar) == 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "fix graphics/labels variable {} is not equal-style variable", pix.xstr);
      pix.xvar = ivar;
    }
    if (pix.ystr) {
      int ivar = input->variable->find(pix.ystr);
      if (ivar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix graphics/labels does not exist", pix.ystr);
      if (input->variable->equalstyle(ivar) == 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "fix graphics/labels variable {} is not equal-style variable", pix.ystr);
      pix.yvar = ivar;
    }
    if (pix.zstr) {
      int ivar = input->variable->find(pix.zstr);
      if (ivar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix graphics/labels does not exist", pix.zstr);
      if (input->variable->equalstyle(ivar) == 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "fix graphics/labels variable {} is not equal-style variable", pix.zstr);
      pix.zvar = ivar;
    }
    if (pix.sstr) {
      int ivar = input->variable->find(pix.sstr);
      if (ivar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix graphics/labels does not exist", pix.sstr);
      if (input->variable->equalstyle(ivar) == 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "fix graphics/labels variable {} is not equal-style variable", pix.sstr);
      pix.svar = ivar;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGraphicsLabels::setup(int vflag)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGraphicsLabels::end_of_step()
{
  numobjs = pixmaps.size();
  if (numobjs == 0) return;

  if (varflag) modify->clearstep_compute();

  if (comm->me == 0) {

    memory->destroy(imgobjs);
    memory->destroy(imgparms);
    memory->create(imgobjs, numobjs, "fix_graphics_labels:imgobjs");
    memory->create(imgparms, numobjs, 11, "fix_graphics_labels:imgparms");

    int n = 0;
    for (auto &pix : pixmaps) {
      if (pix.xstr) pix.pos[0] = input->variable->compute_equal(pix.xvar);
      if (pix.ystr) pix.pos[1] = input->variable->compute_equal(pix.yvar);
      if (pix.zstr) pix.pos[2] = input->variable->compute_equal(pix.zvar);
      if (pix.sstr) pix.scale = input->variable->compute_equal(pix.svar);

      imgobjs[n] = DumpImage::PIXMAP;
      imgparms[n][0] = 1;
      imgparms[n][1] = pix.pos[0];
      imgparms[n][2] = pix.pos[1];
      imgparms[n][3] = pix.pos[2];
      imgparms[n][4] = pix.width;
      imgparms[n][5] = pix.height;
      imgparms[n][6] = ubuf((int64_t) pix.pixmap).d;
      imgparms[n][7] = pix.transcolor[0];
      imgparms[n][8] = pix.transcolor[1];
      imgparms[n][9] = pix.transcolor[2];
      imgparms[n][10] = pix.scale;
      ++n;
    }
  }

  if (varflag) modify->addstep_compute((update->ntimestep / nevery) * nevery + nevery);
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image
------------------------------------------------------------------------- */

int FixGraphicsLabels::image(int *&objs, double **&parms)
{
  if (comm->me == 0) {
    objs = imgobjs;
    parms = imgparms;
    return numobjs;
  }
  return 0;
}
