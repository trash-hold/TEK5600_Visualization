#include "helpers.h"

#include <algorithm>

using namespace H5;

void readH5File(const std::string& filename, RawData* data){
    try {
        // Open file in read-only mode
        H5File file(filename, H5F_ACC_RDONLY);

        // There are two datasets in both files:
        DataSet x_dataset = file.openDataSet("/Velocity/X-comp");
        DataSet y_dataset = file.openDataSet("/Velocity/Y-comp");

        // Getting size of the datasets (assuming both datasets have the same dimensions)
        hsize_t dims[2];
        DataSpace space = x_dataset.getSpace();
        space.getSimpleExtentDims(dims);

        // The read data will be stored as 1D vector -- it's faster
        data->n_rows = dims[0];
        data->n_cols = dims[1];
        data->x.resize(dims[0] * dims[1]);
        data->y.resize(dims[0] * dims[1]);

        x_dataset.read(data->x.data(), PredType::NATIVE_FLOAT);
        y_dataset.read(data->y.data(), PredType::NATIVE_FLOAT);

        std::cout << dims[0] << " x " << dims[1] << "\n";

        file.close();
    }
    catch (H5::Exception& e) {
        std::cerr << "Error reading file: " << e.getCDetailMsg() << std::endl;
    }
}

float bi_interpolate(float v00, float v01, float v10, float v11, float x_frac, float y_frac) {
    // The points are interprated as follows:
    // (v10)---------(v11)
    //    |             |
    //    |        P    |
    //    |         .   |     _
    //    |             |     |
    //    |             |     |  y_frac
    // (v00)---------(v01)    -
    //
    //    |---------|  x_frac
    // 
    // And the billinear interpolation is done according to the formula found in the paper. 
    
    float v0 = v00*(1-x_frac) + v01*x_frac;
    float v1 = v10*(1-x_frac) + v11*x_frac;

    return v0*(1-y_frac) + v1*y_frac;
}

Vec2 bi_interpolate(const RawData* data, const Vec2& pos)
{
    const float eps = 1e-4f;

    // For numerical stability we add -eps to the upper bound of the clamping, to avoid the case when we are exactly on the edge and try to access out of bounds data.
    float x = std::clamp(pos.x, 0.0f, static_cast<float>(data->n_cols - 2) - eps);
    float y = std::clamp(pos.y, 0.0f, static_cast<float>(data->n_rows - 2) - eps);

    size_t col = static_cast<size_t>(x);
    size_t row = static_cast<size_t>(y);

    float x_frac = x - static_cast<float>(col);
    float y_frac = y - static_cast<float>(row);

    Vec2 v00 = data->at(row, col);
    Vec2 v01 = data->at(row, col + 1);
    Vec2 v10 = data->at(row + 1, col);
    Vec2 v11 = data->at(row + 1, col + 1);

    float inter_x = bi_interpolate(v00.x, v01.x, v10.x, v11.x, x_frac, y_frac);
    float inter_y = bi_interpolate(v00.y, v01.y, v10.y, v11.y, x_frac, y_frac);

    return Vec2(inter_x, inter_y);
}

Vec2 bi_interpolate(const PixelPlane* data, const Vec2& pos)
{
    const float eps = 1e-4f;

    // For numerical stability we add -eps to the upper bound of the clamping, to avoid the case when we are exactly on the edge and try to access out of bounds data.
    float x = std::clamp(pos.x, 0.0f, static_cast<float>(data->n_cols - 2) - eps);
    float y = std::clamp(pos.y, 0.0f, static_cast<float>(data->n_rows - 2) - eps);

    size_t col = static_cast<size_t>(x);
    size_t row = static_cast<size_t>(y);

    float x_frac = x - static_cast<float>(col);
    float y_frac = y - static_cast<float>(row);

    Vec2 v00 = data->at(row, col);
    Vec2 v01 = data->at(row, col + 1);
    Vec2 v10 = data->at(row + 1, col);
    Vec2 v11 = data->at(row + 1, col + 1);

    float inter_x = bi_interpolate(v00.x, v01.x, v10.x, v11.x, x_frac, y_frac);
    float inter_y = bi_interpolate(v00.y, v01.y, v10.y, v11.y, x_frac, y_frac);

    return Vec2(inter_x, inter_y);
}