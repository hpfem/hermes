// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "global.h"
#include "quad_all.h"
#include "precalc.h"
#include "mesh.h"
namespace Hermes
{
  namespace Hermes2D
  {
    PrecalcShapeset::PrecalcShapeset(Shapeset* shapeset) : Function<double>()
    {
      if (shapeset == nullptr)
        throw Exceptions::NullException(0);
      this->shapeset = shapeset;
      num_components = shapeset->get_num_components();
      update_max_index();
      set_quad_2d(&g_quad_2d_std);
    }

    void PrecalcShapeset::update_max_index()
    {
      max_index[0] = shapeset->get_max_index(HERMES_MODE_TRIANGLE);
      max_index[1] = shapeset->get_max_index(HERMES_MODE_QUAD);
    }

    void PrecalcShapeset::set_quad_2d(Quad2D* quad_2d)
    {
      Function<double>::set_quad_2d(quad_2d);
    }

    void PrecalcShapeset::set_active_shape(int index)
    {
      // Update the Node table.
      this->invalidate_values();

      this->index = index;
      this->order = std::max(H2D_GET_H_ORDER(shapeset->get_order(index, element->get_mode())), H2D_GET_V_ORDER(shapeset->get_order(index, element->get_mode())));
    }

    void PrecalcShapeset::precalculate(unsigned short order, unsigned short mask)
    {
      Function<double>::precalculate(order, mask);

      unsigned char np = this->quads[cur_quad]->get_num_points(order, this->element->get_mode());
      double3* pt = this->quads[cur_quad]->get_points(order, this->element->get_mode());

      unsigned short j, k;

      ElementMode2D mode = element->get_mode();

      // Correction of points for sub-element mappings.
      if (this->sub_idx != 0)
      {
        for (short i = 0; i < np; i++)
        {
          ref_points[i][0] = ctm->m[0] * pt[i][0] + ctm->t[0];
          ref_points[i][1] = ctm->m[1] * pt[i][1] + ctm->t[1];
        }

        for (j = 0; j < num_components; j++)
          for (k = 0; k < H2D_NUM_FUNCTION_VALUES; k++)
            if (mask & idx2mask[k][j])
              for (short i = 0; i < np; i++)
                this->values[j][k][i] = shapeset->get_value(k, index, ref_points[i][0], ref_points[i][1], j, mode);
      }
      else
      {
        if (this->num_components == 1)
        {
          if (mode == HERMES_MODE_TRIANGLE)
          {
            if (mask & idx2mask[0][0])
              for (short i = 0; i < np; i++)
                this->values[0][0][i] = shapeset->get_fn_value_0_tri(index, pt[i][0], pt[i][1]);
            if (mask & idx2mask[1][0])
              for (short i = 0; i < np; i++)
                this->values[0][1][i] = shapeset->get_dx_value_0_tri(index, pt[i][0], pt[i][1]);
            if (mask & idx2mask[2][0])
              for (short i = 0; i < np; i++)
                this->values[0][2][i] = shapeset->get_dy_value_0_tri(index, pt[i][0], pt[i][1]);
          }
          else
          {
            if (mask & idx2mask[0][0])
              for (short i = 0; i < np; i++)
                this->values[0][0][i] = shapeset->get_fn_value_0_quad(index, pt[i][0], pt[i][1]);
            if (mask & idx2mask[1][0])
              for (short i = 0; i < np; i++)
                this->values[0][1][i] = shapeset->get_dx_value_0_quad(index, pt[i][0], pt[i][1]);
            if (mask & idx2mask[2][0])
              for (short i = 0; i < np; i++)
                this->values[0][2][i] = shapeset->get_dy_value_0_quad(index, pt[i][0], pt[i][1]);
          }
        }
        else
        {
          for (j = 0; j < num_components; j++)
            for (k = 0; k < H2D_NUM_FUNCTION_VALUES; k++)
              if (mask & idx2mask[k][j])
                for (short i = 0; i < np; i++)
                  this->values[j][k][i] = shapeset->get_value(k, index, pt[i][0], pt[i][1], j, mode);
        }
      }
    }

    void PrecalcShapeset::free()
    {
    }

    PrecalcShapeset::~PrecalcShapeset()
    {
      free();
    }

    int PrecalcShapeset::get_active_shape() const
    {
      return index;
    };

    Shapeset* PrecalcShapeset::get_shapeset() const
    {
      return shapeset;
    }

    SpaceType PrecalcShapeset::get_space_type() const
    {
      return shapeset->get_space_type();
    }

    unsigned short PrecalcShapeset::get_edge_fn_order(int edge)
    {
      return H2D_MAKE_EDGE_ORDER(element->get_mode(), edge, shapeset->get_order(index, element->get_mode()));
    }

    std::vector<PrecalcShapesetAssemblingStorage*> PrecalcShapesetAssembling::tables;
    std::vector<PrecalcShapesetAssembling> temp;

    PrecalcShapesetAssembling::PrecalcShapesetAssembling(Shapeset* shapeset) : PrecalcShapeset(shapeset), storage(nullptr)
    {
      for (unsigned short i = 0; i < tables.size(); i++)
      {
        if (tables[i]->shapeset_id == shapeset->get_id())
        {
          this->storage = tables[i];
          this->storage->ref_count++;
          break;
        }
      }
      if (!this->storage)
      {
#pragma omp critical (pss_table_creation)
        {
          this->storage = new PrecalcShapesetAssemblingStorage(this->shapeset);
          this->storage->ref_count++;
          tables.push_back(storage);
          temp.push_back(PrecalcShapesetAssembling(shapeset));
        }
      }
    }

    PrecalcShapesetAssembling::PrecalcShapesetAssembling(const PrecalcShapesetAssembling& other) : PrecalcShapeset(other.shapeset)
    {
      this->storage = other.storage;
      this->storage->ref_count++;
    }

    PrecalcShapesetAssembling::~PrecalcShapesetAssembling()
    {
      if (this->storage && this->storage->ref_count > 0)
      {
#pragma omp atomic
        this->storage->ref_count--;
        if (this->storage && this->storage->ref_count == 0)
        {
#pragma omp critical
          {
            if (this->storage && this->storage->ref_count == 0)
            {
              delete this->storage;
              this->storage = nullptr;
            }
          }
        }
      }
    }

    const double* PrecalcShapesetAssembling::get_fn_values(int component) const
    {
      unsigned char mode = this->element->get_mode();
      if (this->index >= 0 && this->get_quad_2d()->get_id() == 1 && this->sub_idx == 0 && this->num_components == 1 && this->storage->PrecalculatedInfo[this->element->get_mode()][0][this->order][this->index])
        return this->storage->PrecalculatedValues[mode][0][this->order][this->index];
      assert(this->values_valid);
      return &values[component][0][0];
    }

    const double* PrecalcShapesetAssembling::get_dx_values(int component) const
    {
      unsigned char mode = this->element->get_mode();
      if (this->index >= 0 && this->get_quad_2d()->get_id() == 1 && this->sub_idx == 0 && this->num_components == 1 && this->storage->PrecalculatedInfo[this->element->get_mode()][0][this->order][this->index])
        return this->storage->PrecalculatedValues[mode][1][this->order][this->index];
      assert(this->values_valid);
      return &values[component][1][0];
    }

    const double* PrecalcShapesetAssembling::get_dy_values(int component) const
    {
      unsigned char mode = this->element->get_mode();
      if (this->index >= 0 && this->get_quad_2d()->get_id() == 1 && this->sub_idx == 0 && this->num_components == 1 && this->storage->PrecalculatedInfo[this->element->get_mode()][0][this->order][this->index])
        return this->storage->PrecalculatedValues[mode][2][this->order][this->index];
      assert(this->values_valid);
      return &values[component][2][0];
    }

#ifdef H2D_USE_SECOND_DERIVATIVES
    const double* PrecalcShapesetAssembling::get_dxx_values(int component) const
    {
      assert(this->values_valid);
      return &values[component][3][0];
    }

    const double* PrecalcShapesetAssembling::get_dyy_values(int component) const
    {
      assert(this->values_valid);
      return &values[component][4][0];
    }

    const double* PrecalcShapesetAssembling::get_dxy_values(int component) const
    {
      assert(this->values_valid);
      return &values[component][5][0];
    }
#endif

    const double* PrecalcShapesetAssembling::get_values(int component, unsigned short item) const
    {
      if (item == 0)
        return this->get_fn_values(component);
      else if (item == 1)
        return this->get_dx_values(component);
      else if (item == 2)
        return this->get_dy_values(component);
      else
        return Function<double>::get_values(component, item);
    }

    void PrecalcShapesetAssembling::precalculate(unsigned short order, unsigned short mask)
    {
      if (this->index >= 0 && this->get_quad_2d()->get_id() == 1 && this->sub_idx == 0 && this->num_components == 1 && this->storage->PrecalculatedInfo[this->element->get_mode()][0][order][index])
        return;
      else
      {
        Function<double>::precalculate(order, mask);

        unsigned char np = this->quads[cur_quad]->get_num_points(order, this->element->get_mode());
        double3* pt = this->quads[cur_quad]->get_points(order, this->element->get_mode());

        unsigned short j, k;

        ElementMode2D mode = element->get_mode();

        // Correction of points for sub-element mappings.
        if (this->sub_idx != 0)
        {
          for (short i = 0; i < np; i++)
          {
            ref_points[i][0] = ctm->m[0] * pt[i][0] + ctm->t[0];
            ref_points[i][1] = ctm->m[1] * pt[i][1] + ctm->t[1];
          }

          for (j = 0; j < num_components; j++)
            for (k = 0; k < H2D_NUM_FUNCTION_VALUES; k++)
              if (mask & idx2mask[k][j])
                for (short i = 0; i < np; i++)
                  this->values[j][k][i] = shapeset->get_value(k, index, ref_points[i][0], ref_points[i][1], j, mode);
        }
        else
        {
          if (this->num_components == 1)
          {
            if (this->index > 0 && this->get_quad_2d()->get_id() == 1)
            {
#pragma omp critical (precalculatingPSS)
              {
                if (!this->storage->PrecalculatedInfo[mode][0][order][index])
                {
                  double* valuePointer;
                  if (mode == HERMES_MODE_TRIANGLE)
                  {
                    valuePointer = this->storage->PrecalculatedValues[HERMES_MODE_TRIANGLE][0][order][index];
                    for (short i = 0; i < np; i++)
                      valuePointer[i] = shapeset->get_fn_value_0_tri(index, pt[i][0], pt[i][1]);

                    valuePointer = this->storage->PrecalculatedValues[HERMES_MODE_TRIANGLE][1][order][index];
                    for (short i = 0; i < np; i++)
                      valuePointer[i] = shapeset->get_dx_value_0_tri(index, pt[i][0], pt[i][1]);

                    valuePointer = this->storage->PrecalculatedValues[HERMES_MODE_TRIANGLE][2][order][index];
                    for (short i = 0; i < np; i++)
                      valuePointer[i] = shapeset->get_dy_value_0_tri(index, pt[i][0], pt[i][1]);

                    this->storage->PrecalculatedInfo[HERMES_MODE_TRIANGLE][0][order][index] = true;
                  }
                  else
                  {
                    valuePointer = this->storage->PrecalculatedValues[HERMES_MODE_QUAD][0][order][index];
                    for (short i = 0; i < np; i++)
                      valuePointer[i] = shapeset->get_fn_value_0_quad(index, pt[i][0], pt[i][1]);

                    valuePointer = this->storage->PrecalculatedValues[HERMES_MODE_QUAD][1][order][index];
                    for (short i = 0; i < np; i++)
                      valuePointer[i] = shapeset->get_dx_value_0_quad(index, pt[i][0], pt[i][1]);

                    valuePointer = this->storage->PrecalculatedValues[HERMES_MODE_QUAD][2][order][index];
                    for (short i = 0; i < np; i++)
                      valuePointer[i] = shapeset->get_dy_value_0_quad(index, pt[i][0], pt[i][1]);

                    this->storage->PrecalculatedInfo[HERMES_MODE_QUAD][0][order][index] = true;
                  }
                }
              }
            }
            else
            {
              if (mode == HERMES_MODE_TRIANGLE)
              {
                if (mask & idx2mask[0][0])
                  for (short i = 0; i < np; i++)
                    this->values[0][0][i] = shapeset->get_fn_value_0_tri(index, pt[i][0], pt[i][1]);
                if (mask & idx2mask[1][0])
                  for (short i = 0; i < np; i++)
                    this->values[0][1][i] = shapeset->get_dx_value_0_tri(index, pt[i][0], pt[i][1]);
                if (mask & idx2mask[2][0])
                  for (short i = 0; i < np; i++)
                    this->values[0][2][i] = shapeset->get_dy_value_0_tri(index, pt[i][0], pt[i][1]);
              }
              else
              {
                if (mask & idx2mask[0][0])
                  for (short i = 0; i < np; i++)
                    this->values[0][0][i] = shapeset->get_fn_value_0_quad(index, pt[i][0], pt[i][1]);
                if (mask & idx2mask[1][0])
                  for (short i = 0; i < np; i++)
                    this->values[0][1][i] = shapeset->get_dx_value_0_quad(index, pt[i][0], pt[i][1]);
                if (mask & idx2mask[2][0])
                  for (short i = 0; i < np; i++)
                    this->values[0][2][i] = shapeset->get_dy_value_0_quad(index, pt[i][0], pt[i][1]);
              }
            }
          }
          else
          {
            for (j = 0; j < num_components; j++)
              for (k = 0; k < H2D_NUM_FUNCTION_VALUES; k++)
                if (mask & idx2mask[k][j])
                  for (short i = 0; i < np; i++)
                    this->values[j][k][i] = shapeset->get_value(k, index, pt[i][0], pt[i][1], j, mode);
          }
        }
      }
    }

    PrecalcShapesetAssemblingStorage::PrecalcShapesetAssemblingStorage(Shapeset* shapeset) : shapeset_id(shapeset->get_id()), ref_count(0)
    {
      this->max_index[0] = shapeset->get_max_index(HERMES_MODE_TRIANGLE);
      this->max_index[1] = shapeset->get_max_index(HERMES_MODE_QUAD);

      for (int i = 0; i < H2D_NUM_MODES; i++)
      {
        unsigned short g_max, np;
        if (i == HERMES_MODE_TRIANGLE)
        {
          g_max = g_max_tri + 1 + 3 * g_max_tri + 3;
          np = H2D_MAX_INTEGRATION_POINTS_COUNT_TRI;
        }
        else
        {
          g_max = g_max_quad + 1 + 4 * g_max_quad + 4;
          np = H2D_MAX_INTEGRATION_POINTS_COUNT_QUAD;
        }

        unsigned short local_base_size = this->max_index[i] + 1;

        for (int j = 0; j < H2D_NUM_FUNCTION_VALUES; j++)
        {
          this->PrecalculatedValues[i][j] = malloc_with_check<double**>(g_max);
          this->PrecalculatedInfo[i][j] = malloc_with_check<bool*>(g_max);

          for (int k = 0; k < g_max; k++)
          {
            this->PrecalculatedValues[i][j][k] = malloc_with_check<double*>(local_base_size);
            this->PrecalculatedInfo[i][j][k] = calloc_with_check<bool>(local_base_size);

            for (int l = 0; l < local_base_size; l++)
              this->PrecalculatedValues[i][j][k][l] = malloc_with_check<double>(np);
          }
        }
      }
    }

    PrecalcShapesetAssemblingStorage::~PrecalcShapesetAssemblingStorage()
    {
      for (int i = 0; i < H2D_NUM_MODES; i++)
      {
        unsigned short g_max;
        if (i == HERMES_MODE_TRIANGLE)
          g_max = g_max_tri + 1 + 3 * g_max_tri + 3;
        else
          g_max = g_max_quad + 1 + 4 * g_max_quad + 4;

        unsigned short local_base_size = this->max_index[i] + 1;

        for (int j = 0; j < H2D_NUM_FUNCTION_VALUES; j++)
        {
          for (int k = 0; k < g_max; k++)
          {
            for (int l = 0; l < local_base_size; l++)
              free_with_check(this->PrecalculatedValues[i][j][k][l]);
            free_with_check(this->PrecalculatedValues[i][j][k]);
            free_with_check(this->PrecalculatedInfo[i][j][k]);
          }

          free_with_check(this->PrecalculatedValues[i][j]);
          free_with_check(this->PrecalculatedInfo[i][j]);
        }
      }
    }
  }
}