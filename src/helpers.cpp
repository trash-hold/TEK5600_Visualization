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