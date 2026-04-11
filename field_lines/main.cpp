#include <iostream>
#include <vector>
#include <string>
#include <H5Cpp.h>

const std::string isabel_file = "./data/isabel_2d.h5";
const std::string metsim_file = "./data/metsim_2d.h5";

using namespace H5;

struct rawData{
    std::vector<float> x;
    std::vector<float> y;
};


void readH5File(const std::string& filename, rawData* data){
    try {
        // Open file in read-only mode
        H5File file(filename, H5F_ACC_RDONLY);

        // There are two datasets in both files:
        DataSet x_dataset = file.openDataSet("/Velocity/X-comp");
        DataSet y_dataset = file.openDataSet("/Velocity/Y-comp");

        // Getting size of the datasets -> I assume that the vectors are of equal size, so I only check one of them
        hsize_t dims[2];
        DataSpace space = x_dataset.getSpace();
        space.getSimpleExtentDims(dims);

        // The read data will be stored as 1D vector -- it's faster
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

int main(){
    std::cout << "Starting processing of the files..." << std::endl;

    rawData isabel_data;
    rawData metsim_data;

    readH5File(isabel_file, &isabel_data);

    return 0;
}